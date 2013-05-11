/* 
 * Copyright (C) 2011-2013 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * ... and other members of the Reconstruction Team of David Haussler's 
 * lab (BME Dept. UCSC).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE. 
 */
#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafPairCoverageAPI.h"

struct mafCoverageCount {
    // used to slice the coverage data down to the sequence level
    uint64_t sourceLength; // length of the genome / chromosme / sequence in question
    uint64_t count;
};
mafCoverageCount_t* createMafCoverageCount(void) {
    mafCoverageCount_t *mcct = (mafCoverageCount_t *)st_malloc(sizeof(*mcct));
    mcct->sourceLength = 0; // sequence source length
    mcct->count = 0; // number of aligned positions 
    return mcct;
}
uint64_t mafCoverageCount_getSourceLength(mafCoverageCount_t *mcct) {
    return mcct->sourceLength;
}
uint64_t mafCoverageCount_getCount(mafCoverageCount_t *mcct) {
    return mcct->count;
}
void mafCoverageCount_setSourceLength(mafCoverageCount_t *mcct, uint64_t n) {
    mcct->sourceLength = n;
}
void mafCoverageCount_setCount(mafCoverageCount_t *mcct, uint64_t n) {
    mcct->count = n;
}
bool is_wild(const char *s) {
    // return true if char array ends in *, false otherwise
    if (s[strlen(s) - 1] == '*')
        return true;
    return false;
}
bool searchMatched(mafLine_t *ml, const char *seq) {
    // report false if search did not match, true if it did
    if (maf_mafLine_getType(ml) != 's')
        return false;
    if (is_wild(seq)) {
        // only compare up to the wildcard character
        if (!(strncmp(maf_mafLine_getSpecies(ml), seq, strlen(seq) - 1) == 0)) {
            return false;
        }
    } else {
        if (!(strcmp(maf_mafLine_getSpecies(ml), seq) == 0)) {
            return false;
        }
    }
    return true;
}
void compareLines(mafLine_t *ml1, mafLine_t *ml2, stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions) {
    // look through the sequences position by position and count the number of places where the two 
    // sequences contained aligned residues, i.e. neither position contains a gap character.
    char *s1 = maf_mafLine_getSequence(ml1);
    char *s2 = maf_mafLine_getSequence(ml2);
    assert(s1 != NULL);
    assert(s2 != NULL);
    uint64_t n = strlen(s1);
    assert(n == strlen(s2));
    mafCoverageCount_t *mcct1 = stHash_search(seq1Hash, maf_mafLine_getSpecies(ml1));
    mafCoverageCount_t *mcct2 = stHash_search(seq2Hash, maf_mafLine_getSpecies(ml2));
    if (mcct2 == NULL) {
        printf("hash2 is devoid of %s\n", maf_mafLine_getSpecies(ml2));
    }
    assert(mcct1 != NULL);
    assert(mcct2 != NULL);
    for (uint64_t i = 0; i < n; ++i) {
        if ((s1[i] != '-') && (s2[i] != '-')) {
            ++(*alignedPositions);
            ++(mcct1->count);
            ++(mcct2->count);
        }
    }
}
void wrapDestroyMafLine(void *p) {
    maf_destroyMafLineList((mafLine_t *) p);
}
void checkBlock(mafBlock_t *b, const char *seq1, const char *seq2, 
                stHash *seq1Hash, stHash *seq2Hash, uint64_t *alignedPositions) {
    // read through each line of a mafBlock and if the sequence matches the region
    // we're looking for, report the block.
    mafLine_t *ml1 = maf_mafBlock_getHeadLine(b);
    // do a quick scan for either seq before doing the full n^2 comparison
    // im doing this because i know there are some transitively closed mafs that 
    // contain upwards of tens of millions of rows... :/
    bool has1 = false, has2 = false;
    stList *seq1List = stList_construct3(0, wrapDestroyMafLine);
    stList *seq2List = stList_construct3(0, wrapDestroyMafLine);
    mafCoverageCount_t *mcct1 = NULL, *mcct2 = NULL;
    while (ml1 != NULL) {
        if (searchMatched(ml1, seq1)) {
            has1 = true;
            stList_append(seq1List, maf_copyMafLine(ml1));
            // create an item in the hash for this sequence
            mcct1 = NULL;
            if ((mcct1 = stHash_search(seq1Hash, maf_mafLine_getSpecies(ml1))) == NULL) {
                // new sequence, add to the hash
                mcct1 = createMafCoverageCount();
                mcct1->sourceLength = maf_mafLine_getSourceLength(ml1);
                stHash_insert(seq1Hash, stString_copy(maf_mafLine_getSpecies(ml1)), mcct1);
            } else {
                assert(mcct1->sourceLength == maf_mafLine_getSourceLength(ml1));
            }
        }
        if (searchMatched(ml1, seq2)) {
            has2 = true;
           stList_append(seq2List, maf_copyMafLine(ml1));
            // create an item in the hash for this sequence
            mcct2 = NULL;
            if ((mcct2 = stHash_search(seq2Hash, maf_mafLine_getSpecies(ml1))) == NULL) {
                // new sequence, add to the hash
                mcct2 = createMafCoverageCount();
                mcct2->sourceLength = maf_mafLine_getSourceLength(ml1);
                stHash_insert(seq2Hash, stString_copy(maf_mafLine_getSpecies(ml1)), mcct2);
            } else {
                assert(mcct2->sourceLength == maf_mafLine_getSourceLength(ml1));
            }
        }
        ml1 = maf_mafLine_getNext(ml1);
    }
    if (!(has1 && has2)) {
        stList_destruct(seq1List);
        stList_destruct(seq2List);
        return;
    }
    // perform the full n^2 scan
    ml1 = NULL;
    mafLine_t *ml2 = NULL;
    stListIterator *slit1 = stList_getIterator(seq1List);
    stListIterator *slit2 = NULL;
    while ((ml1 = stList_getNext(slit1)) != NULL) {
        slit2 = stList_getIterator(seq2List);
        while ((ml2 = stList_getNext(slit2)) != NULL) {
            compareLines(ml1, ml2, seq1Hash, seq2Hash, alignedPositions);
        }
        stList_destructIterator(slit2);
    } 
    stList_destructIterator(slit1);
    stList_destruct(seq1List);
    stList_destruct(seq2List);
}
void processBody(mafFileApi_t *mfa, char *seq1, char *seq2, stHash *seq1Hash, stHash *seq2Hash,
                 uint64_t *alignedPositions) {
    mafBlock_t *thisBlock = NULL;
    *alignedPositions = 0;
    while ((thisBlock = maf_readBlock(mfa)) != NULL) {
        checkBlock(thisBlock, seq1, seq2, seq1Hash, seq2Hash, alignedPositions);
        maf_destroyMafBlockList(thisBlock);
    }
}
