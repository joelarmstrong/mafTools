##############################
# Copyright (C) 2009-2014 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
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
import xml.etree.ElementTree as ET
import unittest
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../lib/')))
import mafToolsTest as mtt

g_speciesTree = "((A,B)C,D)E;"

g_headers = ['''##maf version=1
''']

# (maf 1, maf 2, earlyCoalescences, lateCoalescences, identicalCoalescences)
# All species trees are ((A,B)C,D)E
# More tests are in the "only leaves" category because manually
# figuring out which ancestral coalescences have changed gets hairy
# quickly.
knownValues = [
    # Test 1 -- two identical mafs.
    ("""a score=0.000 tree="((A,B)C,D)E;"
s A 0 1 + 20 A
s B 0 1 + 20 A
s C 0 1 + 20 A
s D 0 1 + 20 A
s E 0 1 + 20 A
""",
     """a tree="((B,A)C,D)E;"
s B 0 1 + 20 A
s A 0 1 + 20 A
# random comment that should be ignored!
s C 0 1 + 20 A
s D 0 1 + 20 A
s E 0 1 + 20 A
""",
     0, 0, 10)]

# (maf 1, maf 2, earlyCoalescences, lateCoalescences, identicalCoalescences)
# All species trees are ((A,B)C,D)E
knownValuesOnlyLeaves = [
    # Test 1 -- two identical mafs.
    ("""a score=0.000 tree="((A,B)C,D)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s D 0 4 + 20 ACTG
""",
     """a tree="((A,B)C,D)E;"
# random comment that should be ignored!
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s D 0 4 + 20 ACTG
""",
     0, 0, 12),
    # Test 2 -- maf 2 completely earlier than maf 1.
    ("""a score=0.000 tree="((A)C,(B)C)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
""",
     """a score=0.000 tree="((A,B)C)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
""",
     4, 0, 0),
    # Test 3 -- maf 2 completely later than maf 1.
    ("""a score=0.000 tree="((A,B)C)E;"
s A 4 8 + 20 ACTG
s B 4 8 + 20 ACTG
""",
     """a score=0.000 tree="((A)C,(B)C)E;"
s A 4 8 + 20 ACTG
s B 4 8 + 20 ACTG
""",
     0, 4, 0),
    # Test 4 -- blocks from tests 2 + 3 concatenated together (sanity
    # check that blocks are considered independently)
    ("""a score=0.000 tree="((A)C,(B)C)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG

a score=0.000 tree="((A,B)C)E;"
s A 4 8 + 20 ACTG
s B 4 8 + 20 ACTG
""",
     """a score=0.000 tree="((A,B)C)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG

a score=0.000 tree="((A)C,(B)C)E;"
s A 4 8 + 20 ACTG
s B 4 8 + 20 ACTG
""",
     4, 4, 0),
    # Test 5 -- check handling of gaps
    ("""a score=0.000 tree="((A,B)C)E;"
s A 0 2 + 20 A--G
s B 0 4 + 20 ACTG
""",
     """a score=0.000 tree="((A,B)C)E;"
s A 0 1 + 20 A
s B 0 1 + 20 A

a tree="((A)C, (B)C)E;"
s A 1 1 + 20 G
s B 3 1 + 20 G
""",
     0, 1, 1),
    # Test 6 -- Only A and D are comparable, and their induced
    # tree should be the same, despite other parts of the tree
    # changing.
    ("""a score=0.000 tree="((A,B)C, D)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s D 0 4 + 20 ACTG
""",
     """a score=0.000 tree="((A)C, D)E;"
s A 0 4 + 20 ACTG
s D 0 4 + 20 ACTG
""",
     0, 0, 4),
    # Test 7 -- A, B, D are all comparable. A,B relationship changes,
    # but not B,D or A,D.
    ("""a score=0.000 tree="((A,B)C,D)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s D 0 4 + 20 ACTG
""",
     """a score=0.000 tree="((A)C, (B)C, D)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s D 0 4 + 20 ACTG
""",
     0, 4, 8),
    # Test 8 -- duplicates in the same sequence. A,B1 becomes later,
    # B1,B2 becomes earlier, A,B2 stays the same.
    ("""a score=0.000 tree="((A,B)C,(B)C)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s B 4 4 + 20 ACTG
""",
     """a score=0.000 tree="((A)C, (B,B)C)E;"
s A 0 4 + 20 ACTG
s B 0 4 + 20 ACTG
s B 4 4 + 20 ACTG
""",
     4, 4, 4)
]

def getAggregateCoalescenceResult(filename, attrib):
    tree = ET.parse(filename)
    return int(tree.find("aggregateCoalescenceResults").attrib[attrib])

class CoalescencesTest(unittest.TestCase):
    def test_coalescences_knownValues_onlyLeaves(self):
        """Test phylogeny comparison using coalescences, on MAFs that only have alignments for the leaves of the block trees."""
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir("coalescencesOnlyLeaves"))
        for i, (maf1, maf2, earlyCoalescences, lateCoalescences, identicalCoalescences) in enumerate(knownValuesOnlyLeaves):
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPhyloComparator')),
                   '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--numSamples=10000', '--logLevel=critical',
                   "--speciesTree=%s" % g_speciesTree,
                   "--onlyLeaves"
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)

            outputFile = os.path.abspath(os.path.join(tmpDir, 'output.xml'))
            self.assertEqual(getAggregateCoalescenceResult(outputFile, 'earlyCoalescences'), earlyCoalescences, "early coalescences don't match on test %d" % (i + 1))
            self.assertEqual(getAggregateCoalescenceResult(outputFile, 'lateCoalescences'), lateCoalescences, "late coalescences don't match on test %d" % (i + 1))
            self.assertEqual(getAggregateCoalescenceResult(outputFile, 'identicalCoalescences'), identicalCoalescences, "identical coalescences don't match on test %d" % (i + 1))
        mtt.removeDir(tmpDir)

    def test_coalescences_knownValues_withAncestors(self):
        """Test phylogeny comparison using coalescences, on MAFs that only have alignments for ancestors and leaves."""
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir("coalescencesWithAncestors"))
        for i, (maf1, maf2, earlyCoalescences, lateCoalescences, identicalCoalescences) in enumerate(knownValues):
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = [os.path.abspath(os.path.join(parent, 'test', 'mafPhyloComparator')),
                   '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--numSamples=10000', '--logLevel=critical',
                   "--speciesTree=%s" % g_speciesTree,
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)

            outputFile = os.path.abspath(os.path.join(tmpDir, 'output.xml'))
            self.assertEqual(getAggregateCoalescenceResult(outputFile, 'earlyCoalescences'), earlyCoalescences, "early coalescences don't match on test %d" % (i + 1))
            self.assertEqual(getAggregateCoalescenceResult(outputFile, 'lateCoalescences'), lateCoalescences, "late coalescences don't match on test %d" % (i + 1))
            self.assertEqual(getAggregateCoalescenceResult(outputFile, 'identicalCoalescences'), identicalCoalescences, "identical coalescences don't match on test %d" % (i + 1))
        mtt.removeDir(tmpDir)

    def test_coalescences_memoryTest(self):
        """Test that valgrind doesn't catch any memory errors when running mafPhyloComparator in coalescence mode."""
        valgrind = mtt.which('valgrind')
        if valgrind is None:
            return
        mtt.makeTempDirParent()
        tmpDir = os.path.abspath(mtt.makeTempDir("coalescencesMem"))

        # only leaves

        for i, (maf1, maf2, earlyCoalescences, lateCoalescences, identicalCoalescences) in enumerate(knownValuesOnlyLeaves):
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafPhyloComparator')),
                   '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--numSamples=10000', '--logLevel=critical',
                   "--speciesTree=%s" % g_speciesTree,
                   "--onlyLeaves"
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)

        # with ancestors
        for i, (maf1, maf2, earlyCoalescences, lateCoalescences, identicalCoalescences) in enumerate(knownValues):
            testMaf1 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf1.maf')), 
                                    maf1, g_headers)
            testMaf2 = mtt.testFile(os.path.abspath(os.path.join(tmpDir, 'maf2.maf')), 
                                    maf2, g_headers)
            parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cmd = mtt.genericValgrind(tmpDir)
            cmd += [os.path.abspath(os.path.join(parent, 'test', 'mafPhyloComparator')),
                   '--mafFile1', os.path.abspath(os.path.join(tmpDir, 'maf1.maf')),
                   '--mafFile2', os.path.abspath(os.path.join(tmpDir, 'maf2.maf')),
                   '--out', os.path.abspath(os.path.join(tmpDir, 'output.xml')),
                   '--numSamples=10000', '--logLevel=critical',
                   "--speciesTree=%s" % g_speciesTree,
                   ]
            mtt.recordCommands([cmd], tmpDir)
            mtt.runCommandsS([cmd], tmpDir)
            passed = mtt.noMemoryErrors(os.path.join(tmpDir, 'valgrind.xml'))
            self.assertTrue(passed)
        mtt.removeDir(tmpDir)
        

if __name__ == '__main__':
    unittest.main()
