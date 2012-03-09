#!/usr/bin/env python2.6
"""
mafValidator
10 Oct 2011
dent earl, dearl (a) soe ucsc edu

Script to validate Multpile Alignment Format (maf) 
files.

"""
##############################
# Copyright (C) 2009-2012 by
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
#
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC).
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
from optparse import OptionParser
import os
import re
import sys

class ValidatorError(Exception): pass
class SourceLengthError(ValidatorError): pass
class SpeciesFieldError(ValidatorError): pass
class MissingAlignmentBlockLineError(ValidatorError): pass
class AlignmentBlockLineKeyValuePairError(ValidatorError): pass
class AlignmentLengthError(ValidatorError): pass
class FieldNumberError(ValidatorError): pass
class FooterError(ValidatorError): pass
class StrandCharacterError(ValidatorError): pass
class StartFieldError(ValidatorError): pass
class SourceSizeFieldError(ValidatorError): pass
class OutOfRangeError(ValidatorError): pass
class HeaderError(ValidatorError): pass
class ILineFormatError(ValidatorError): pass

def initOptions(parser):
   parser.add_option('--maf', dest = 'filename', 
                     help = 'path to maf file to validate.')
   parser.add_option('--testChromNames', dest = 'testChromNames', 
                     action = 'store_true',
                     default = False,
                     help = ('Test that all species fields contain chrom name, i.e.: '
                             's hg19.chr1 ... default = %default'))
   
def checkOptions(options, args, parser):
   if options.filename is None:
      parser.error('specify --maf')
   if not os.path.exists(options.filename):
      parser.error('--maf %s does not exist.' % options.filename)

def validateMaf(filename, testChromNames = False):
   """ returns true on valid maf file
   """ 
   nameRegex = r'(.+?)\.(chr.+)'
   namePat = re.compile(nameRegex)
   f = open(filename, 'r')
   header = f.next()
   validateHeader(header, filename)
   sources = {}
   prevLineWasAlignmentBlock = False
   alignmentFieldLength = None
   prevline = ''
   for lineno, line in enumerate(f, 2):
      line = line.strip()
      
      if line.startswith('#'):
         prevLineWasAlignmentBlock = False
         alignmentFieldLength = None
         prevline = line
         continue
      elif line.startswith('a'):
         validateAlignmentLine(lineno, line, filename)
         prevLineWasAlignmentBlock = True
         alignmentFieldLength = None
      elif line.startswith('s'):
         if not prevLineWasAlignmentBlock:
            raise MissingAlignmentBlockLineError('maf %s has a sequence line that was not preceded '
                                                 'by an alignment line on line number %d: %s' 
                                                 % (filename, lineno, line))
         species, chrom, length, alFieldLen = validateSeqLine(namePat, testChromNames, 
                                                              lineno, line, filename)
         if alignmentFieldLength is None:
            alignmentFieldLength = alFieldLen
         else:
            if alignmentFieldLength != alFieldLen:
               raise AlignmentLengthError('maf %s has a sequence line with an alignment field of different '
                                          'length than the other sequences in the block on line number %d: %s'
                                          % (filename, lineno, line))
         if (species, chrom)  in sources:
            if sources[(species, chrom)] != length:
               raise SourceLengthError('maf %s has different source lengths for '
                                       'species %s lineno %d: %s'
                                       % (filename, species, lineno, line))
         else:
            sources[(species, chrom)] = length
      elif line.startswith('i'):
         if not prevLineWasAlignmentBlock:
            raise MissingAlignmentBlockLineError('maf %s has a sequence line that was not preceded '
                                                 'by an alignment line on line number %d: %s' 
                                                 % (filename, lineno, line))
         validateILine(lineno, line, prevline, filename)
      elif line.startswith('e'):
         if not prevLineWasAlignmentBlock:
            raise MissingAlignmentBlockLineError('maf %s has a sequence line that was not preceded '
                                                 'by an alignment line on line number %d: %s' 
                                                 % (filename, lineno, line))
      elif line.startswith('q'):
         if not prevLineWasAlignmentBlock:
            raise MissingAlignmentBlockLineError('maf %s has a sequence line that was not preceded '
                                                 'by an alignment line on line number %d: %s' 
                                                 % (filename, lineno, line))
      elif line == '':
         prevLineWasAlignmentBlock = False
         alignmentFieldLength = None
      prevline = line

   if line != '':
      raise FooterError('maf %s has a bad footer, should end with two new lines.' 
                        % filename)
   return True

def validateILine(lineno, line, prevline, filename):
   """ Checks all lines that start with 'i' and raises an exepction if
   the line is malformed.
   i lines are made up of five fields after the 'i':
   From http://genome.ucsc.edu/FAQ/FAQformat#format5 :
   src -- The name of the source sequence for the alignment. 
          Should be the same as the 's' line immediately above this line.
   leftStatus -- A character that specifies the relationship between the 
                 sequence in this block and the sequence that appears in 
                 the previous block.
   leftCount -- Usually the number of bases in the aligning species between 
                the start of this alignment and the end of the previous one.
   rightStatus -- A character that specifies the relationship between the 
                  sequence in this block and the sequence that appears in 
                  the subsequent block.
   rightCount -- Usually the number of bases in the aligning species between 
                 the end of this alignment and the start of the next one.
   The status characters can be one of the following values:
   C -- the sequence before or after is contiguous with this block.
   I -- there are bases between the bases in this block and the one before or after it.
   N -- this is the first sequence from this src chrom or scaffold.
   n -- this is the first sequence from this src chrom or scaffold but it is 
        bridged by another alignment from a different chrom or scaffold.
   M -- there is missing data before or after this block (Ns in the sequence).
   T -- the sequence in this block has been used before in a previous block (likely a tandem duplication)
   """
   d = line.split()
   p = prevline.split()
   if len(d) != 6:
      raise ILineFormatError('maf %s contains an "i" line that has too many fields on line number %d: '
                             '%s' % (filename, lineno, line))
   if p[0] != 's':
      raise ILineFormatError('maf %s contains an "i" line that does not follow an "s" line on line number %d: '
                             '%s' % (filename, lineno, line))
   for i in [3, 5]:
      try:
         n = int(d[i])
      except ValueError:
         raise ILineFormatError('maf %s contains an "i" line that has non integer Count "%s" on line number %d: '
                                '%s' % (filename, d[i], lineno, line))
      if int(d[i]) < 0:
         raise ILineFormatError('maf %s contains an "i" line that has negative Count "%s" on line number %d: '
                                '%s' % (filename, d[i], lineno, line))
   for i in [2, 4]:
      if d[i] not in ['C', 'I', 'N', 'n', 'M', 'T']:
         raise ILineFormatError('maf %s contains an "i" line with an invalid Status "%s" on line number %d: '
                                '%s' % (filename, d[i], lineno, line))
      if d[i] == 'I' and int(d[i + 1]) < 1:
         raise ILineFormatError('maf %s contains an "i" line with an invalid Count "%s" on line number %d: '
                                '%s' % (filename, d[i + 1], lineno, line))
   if p[1] != d[1]:
      raise ILineFormatError('maf %s contains an "i" line with a different src value "%s" than on the previous '
                             '"s" line "%s" on line number %d: '
                             '%s' % (filename, d[1], p[1], lineno, line))

def validateAlignmentLine(lineno, line, filename):
   """ Checks all lines that start with 'a' and raises an exception if 
   the line is malformed.
   """
   d = line.split()
   for i in xrange(1, len(d)):
      if len(d[i].split('=')) != 2:
         raise AlignmentBlockLineKeyValuePairError('maf %s has an alignment line that does not contain '
                                                   'good key-value pairs on line number %d: %s' 
                                                   % (filename, lineno, line))

def validateSeqLine(namePat, testChromNames, lineno, line, filename):
   data = line.split()
   if len(data) != 7:
      raise FieldNumberError('maf %s has incorrect number of fields on line number %d: %s' 
                             % (filename, lineno, line))
   if data[4] not in ('-', '+'):
      raise StrandCharacterError('maf %s has unexpected character in strand field "%s" lineno %d: %s' 
                                 % (filename, data[4], lineno, line))
   if int(data[3]) != len(data[6].replace('-', '')):
      raise AlignmentLengthError('maf %s has incorrect seq len (should be %d) or alignment field lineno %d: %s'
                                 % (filename, len(data[6].replace('-', '')), lineno, line))
   if int(data[2]) < 0:
      raise StartFieldError('maf %s has bad start field lineno %d: %s'
                           % (filename, lineno, line))
   if int(data[5]) < 0:
      raise SourceSizeFieldError('maf %s has bad srcSize field lineno %d: %s'
                                 % (filename, lineno, line))
   if int(data[2]) + int(data[3]) > int(data[5]):
      raise OutOfRangeError('maf %s out of range sequence lineno %d: %s'
                            % (filename, lineno, line))
   if testChromNames:
      m = re.match(namePat, data[1])
      if m is None:
         raise SpeciesFieldError('maf %s has name (src) field without ".chr" suffix: "%s" lineno %d: %s' 
                                 % (filename, data[1], lineno, line))
      return m.group(1), m.group(2), data[5]
   return data[1], None, data[5], len(data[6])

def validateHeader(header, filename):
   """ tests the first line of the maf file to make sure it is valid
   """
   if not header.startswith('##'):
      raise HeaderError('maf %s has bad header, fails to start with `##\': %s' 
                           % (filename, header))
   data = header.split()
   version = False
   for d in data[1:]:
      if d.startswith('=') or d == '=' or d.endswith('='):
         raise HeaderError('maf %s has bad header, there may be '
                              'no whitespace surrounding "=": %s' 
                              % (filename, header))
      try:
         k,v = d.split('=')
      except ValueError:
         raise HeaderError('maf %s has bad header, there may be '
                              'no whitespace surrounding "=": %s' 
                              % (filename, header))
      if k == 'version':
         version = True
   if not version:
      raise HeaderError('maf %s has bad header, no version information: %s' 
                           % (filename, header))
         
def main():
   parser = OptionParser()
   initOptions(parser)
   options, args = parser.parse_args()
   checkOptions(options, args, parser)
   
   validateMaf(options.filename, options.testChromNames)

if __name__ == '__main__':
   main()
