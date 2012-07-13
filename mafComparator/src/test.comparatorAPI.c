/* 
 * Copyright (C) 2012 by 
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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "CuTest.h"
#include "common.h"
#include "sonLib.h"
#include "comparatorAPI.h"

static void printMat(uint64_t **mat, uint64_t n) {
    printf("printMat(mat, %" PRIi64 ")\n", n);
    for (uint64_t r = 0; r < n; ++r) {
        for (uint64_t c = 0; c < n; ++c) {
            if (r == c) {
                printf("  - ");
            } else if (c < r) {
                printf("    ");
            } else {
                printf(" %2" PRIi64 " ", mat[r][c]);
            }
        }
        printf("\n");
    }
    printf("\n");
}
static void test_mappingMatrixToArray_0(CuTest *testCase) {
    uint64_t **mat = NULL;
    uint64_t index, result;
    for (uint64_t n = 2; n < 11; ++n) {
        // build out mat
        mat = (uint64_t **) st_malloc(sizeof(*mat) * n);
        for (uint64_t r = 0; r < n; ++r) {
            mat[r] = st_malloc(sizeof(uint64_t *) * n);
        }
        index = 0;
        for (uint64_t r = 0; r < n; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                mat[r][c] = index++;
            }
        }
        // printMat(mat, n);
        // test
        for (uint64_t r = 0; r < n - 1; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                pairIndicesToArrayIndex(r, c, n, &result);
                CuAssertTrue(testCase, mat[r][c] == result);
            }
        }
        // clean up
        for (uint64_t r = 0; r < n; ++r) {
            free(mat[r]);
        }
        free(mat);
    }
}
static void test_mappingArrayToMatrix_0(CuTest *testCase) {
    uint64_t **mat = NULL;
    uint64_t index, resultRow, resultCol;
    uint64_t n;
    for (uint64_t p = 0; p < 12; ++p) {
        n = 2 << p;
        // build out mat
        mat = (uint64_t **) st_malloc(sizeof(*mat) * n);
        for (uint64_t r = 0; r < n; ++r) {
            mat[r] = st_malloc(sizeof(uint64_t *) * n);
        }
        index = 0;
        for (uint64_t r = 0; r < n; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                mat[r][c] = index++;
            }
        }
        // printMat(mat, n);
        // test
        for (uint64_t r = 0; r < n - 1; ++r) {
            for (uint64_t c = r + 1; c < n; ++c) {
                arrayIndexToPairIndices(mat[r][c], n, &resultRow, &resultCol);
                CuAssertTrue(testCase, r == resultRow);
                CuAssertTrue(testCase, c == resultCol);
            }
        }
        // clean up
        for (uint64_t r = 0; r < n; ++r) {
            free(mat[r]);
        }
        free(mat);
    }
}
static void test_pairCounting_0(CuTest *testCase) {
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    uint64_t expected, observed;
    uint32_t lineNumber, numSeqs, numLines;
    stHash *legitPairs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    mafLine_t *ml = NULL;
    mafBlock_t *mb = NULL;
    numSeqs = 4;
    numLines = 5;
    expected = chooseTwo(4) * 13;
    char **input = (char**) st_malloc(sizeof(*input) * numLines);
    input[0] = stString_copy("a score=0.0");
    input[1] = stString_copy("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca");
    input[2] = stString_copy("s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca");
    input[3] = stString_copy("s baboon         249182 13 +   4622798 gcagctgaaaaca");
    input[4] = stString_copy("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA");
    stHash_insert(legitPairs, stString_copy("hg16.chr7"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("panTro1.chr6"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("baboon"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("mm4.chr6"), stString_copy(""));
    mb = maf_newMafBlock();
    lineNumber = 3;
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_setNumberOfSequences(mb, numSeqs);
    maf_mafBlock_setNumberOfLines(mb, numLines);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString(input[0], lineNumber));
    ml = maf_mafBlock_getHeadLine(mb);
    for (uint32_t i = 0; i < maf_mafBlock_getNumberOfLines(mb); ++i, ++lineNumber) {
        maf_mafLine_setNext(ml, maf_newMafLineFromString(input[i], lineNumber));
        ml = maf_mafLine_getNext(ml);
        maf_mafBlock_setTailLine(mb, ml);
    }
    observed = walkBlockCountingPairs(mb, legitPairs, chooseTwoArray);
    CuAssertTrue(testCase, observed == expected);
    // clean up
    for (unsigned i = 0; i < numLines; ++i) { 
        free(input[i]);
    }
    free(input);
    maf_destroyMafBlockList(mb);
    stHash_destruct(legitPairs);

    // new test
    expected = chooseTwo(3) * 13;
    legitPairs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    numSeqs = 4;
    numLines = 5;
    input = (char**) st_malloc(sizeof(*input) * numLines);
    input[0] = stString_copy("a score=0.0");
    input[1] = stString_copy("s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca");
    input[2] = stString_copy("s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca");
    input[3] = stString_copy("s baboon         249182 13 +   4622798 gcagctgaaaaca");
    input[4] = stString_copy("s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA");
    stHash_insert(legitPairs, stString_copy("hg16.chr7"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("panTro1.chr6"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("mm4.chr6"), stString_copy(""));
    mb = maf_newMafBlock();
    lineNumber = 3;
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_setNumberOfSequences(mb, numSeqs);
    maf_mafBlock_setNumberOfLines(mb, numLines);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString(input[0], lineNumber));
    ml = maf_mafBlock_getHeadLine(mb);
    for (uint32_t i = 0; i < maf_mafBlock_getNumberOfLines(mb); ++i, ++lineNumber) {
        maf_mafLine_setNext(ml, maf_newMafLineFromString(input[i], lineNumber));
        ml = maf_mafLine_getNext(ml);
        maf_mafBlock_setTailLine(mb, ml);
    }
    observed = walkBlockCountingPairs(mb, legitPairs, chooseTwoArray);
    CuAssertTrue(testCase, observed == expected);
    // clean up
    for (unsigned i = 0; i < numLines; ++i) { 
        free(input[i]);
    }
    free(input);
    maf_destroyMafBlockList(mb);
    stHash_destruct(legitPairs);
    // new test
    expected = chooseTwo(3) + chooseTwo(4) * 37;
    legitPairs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    numSeqs = 5;
    numLines = 6;
    input = (char**) st_malloc(sizeof(*input) * numLines);
    input[0] = stString_copy("a score=23262.0   ");
    input[1] = stString_copy("s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    input[2] = stString_copy("s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG");
    input[3] = stString_copy("s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG");
    input[4] = stString_copy("s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG");
    input[5] = stString_copy("s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG");
    stHash_insert(legitPairs, stString_copy("hg18.chr7"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("panTro1.chr6"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("baboon"), stString_copy(""));
    stHash_insert(legitPairs, stString_copy("rn3.chr4"), stString_copy(""));
    mb = maf_newMafBlock();
    lineNumber = 3;
    maf_mafBlock_setLineNumber(mb, lineNumber);
    maf_mafBlock_setNumberOfSequences(mb, numSeqs);
    maf_mafBlock_setNumberOfLines(mb, numLines);
    maf_mafBlock_setHeadLine(mb, maf_newMafLineFromString(input[0], lineNumber));
    ml = maf_mafBlock_getHeadLine(mb);
    for (uint32_t i = 0; i < maf_mafBlock_getNumberOfLines(mb); ++i, ++lineNumber) {
        maf_mafLine_setNext(ml, maf_newMafLineFromString(input[i], lineNumber));
        ml = maf_mafLine_getNext(ml);
        maf_mafBlock_setTailLine(mb, ml);
    }
    observed = walkBlockCountingPairs(mb, legitPairs, chooseTwoArray);
    CuAssertTrue(testCase, observed == expected);
    // clean up
    for (unsigned i = 0; i < numLines; ++i) { 
        free(input[i]);
    }
    free(input);
    free(chooseTwoArray);
    maf_destroyMafBlockList(mb);
    stHash_destruct(legitPairs);
}
static char** createRandomColumn(uint32_t n, uint32_t colLength, double gapProb) {
    char **mat = (char**) st_malloc(sizeof(*mat) * n);
    for (uint32_t i = 0; i < n; ++i) {
        mat[i] = st_malloc(sizeof(char) * colLength + 1);
        for (uint32_t j = 0; j < colLength; ++j) {
            if (st_random() <= gapProb) {
                mat[i][j] = '-';
            } else {
                mat[i][j] = 'A';
            }
        }
        mat[i][colLength] = '\0'; // unused
    }
    return mat;
}
static bool* createRandomLegitRow(uint32_t n, double alpha) {
    bool *legitRow = (bool*) st_malloc(sizeof(*legitRow) * n);
    for (uint32_t i = 0; i < n; ++i) {
        if (st_random() <= alpha) {
            legitRow[i] = true;
        } else {
            legitRow[i] = false;
        }
    }
    return legitRow;
}
static double runningAverage(double m, double x, uint32_t i) {
    m += (x - m) / (i + 1);
    return m;
}
static char* randomName(uint32_t n) {
    char *s = (char*) st_malloc(sizeof(*s) * (n + 9));
    strcpy(s, "species_");
    for (uint32_t i = 8; i < n + 8; ++i) {
        if (st_random() < 0.5) {
            // upper case
            s[i] = st_randomInt(65, 91);
        } else {
            // lower case 
            s[i] = st_randomInt(97, 123);
        }
    }
    s[n] = '\0';
    return s;
}
static mafLine_t** createMlArray(uint32_t n) {
    mafLine_t **mlArray = (mafLine_t**) st_malloc(sizeof(*mlArray) * n);
    for (uint32_t i = 0; i < n; ++i) {
        mlArray[i] = maf_newMafLine();
        maf_mafLine_setSpecies(mlArray[i], randomName(20));
    }
    return mlArray;
}
static char** createNameArray(uint32_t n) {
    char **nameArray = (char**) st_malloc(sizeof(*nameArray) * n);
    for (uint32_t i = 0; i < n; ++i) {
        nameArray[i] = randomName(20); // this now needs to be free'd
    }
    return nameArray;
}
static void u32set(uint32_t *a, uint64_t v, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        a[i] = v;
    }
}
static void intset(int *a, int v, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        a[i] = v;
    }
}
static void test_columnSampling_timing_0(CuTest *testCase) {
    // this should not be enabled for the default test.
    (void) (createMlArray);
    char **mat = NULL;
    bool *legitRows = NULL;
    uint32_t n, m;
    uint32_t *positions = NULL;
    int *strandInts = NULL;
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    double timeClever, timeNaive, p;
    time_t t1;
    uint32_t colLength = 2000;
    stSortedSet *pairs = NULL;
    char **nameArray = NULL;
    printf("#Rows        p      n*p clever naive\n");
    for (uint32_t i = 0; i < 9; ++i) {
        pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, 
                                       (void(*)(void *)) aPair_destruct);
        n = 2 << i;
        p = 2.0 / (n * (n - 1));
        mat = createRandomColumn(n, colLength, 0.1);
        legitRows = createRandomLegitRow(n, 0.9);
        m = sumBoolArray(legitRows, n);
        positions = (uint32_t*) st_malloc(sizeof(*positions) * n);
        u32set(positions, 0, n);
        strandInts = (int*) st_malloc(sizeof(*strandInts) * n);
        intset(strandInts, 1, n);
        nameArray = createNameArray(n);
        timeClever = 0.0;
        timeNaive = 0.0;
        t1 = time(NULL);
        for (uint32_t c = 0; c < colLength; ++c) {
            samplePairsFromColumn(0.01, pairs, m, chooseTwoArray, nameArray, positions);
            updatePositions(mat, c, positions, strandInts, n);
        }
        timeClever = difftime(time(NULL), t1);
        free(positions);
        stSortedSet_destruct(pairs);
        positions = (uint32_t*) st_malloc(sizeof(*positions) * n);
        memset(positions, 0, sizeof(*positions) * n);
        pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, 
                                       (void(*)(void *)) aPair_destruct);
        t1 = time(NULL);
        for (uint32_t c = 0; c < colLength; ++c) {
            samplePairsFromColumnNaive(mat, c, legitRows, 0.01, pairs, chooseTwoArray, 
                                       nameArray, positions, n, chooseTwo(n));
            updatePositions(mat, c, positions, strandInts, n);
        }
        timeNaive = difftime(time(NULL), t1);
        printf("%5" PRIu32 " %6.2e %6f %4.0fs %4.0fs\n", n, p, p * n, timeClever, timeNaive);
        // clean up
        for (uint32_t j = 0; j < n; ++j) {
            free(mat[j]);
            free(nameArray[j]); // we ONLY do this in this test example, not in production code.
        }
        free(mat);
        free(nameArray);
        free(legitRows);
    }
    free(chooseTwoArray);
    CuAssertTrue(testCase, true);
}
static void test_chooseTwoValues_0(CuTest *testCase) {
    CuAssertTrue(testCase, chooseTwo(0) == 0);
    CuAssertTrue(testCase, chooseTwo(1) == 0);
    CuAssertTrue(testCase, chooseTwo(2) == 1);
    CuAssertTrue(testCase, chooseTwo(3) == 3);
    CuAssertTrue(testCase, chooseTwo(4) == 6);
    CuAssertTrue(testCase, chooseTwo(5) == 10);
    CuAssertTrue(testCase, chooseTwo(48) == 1128);
    CuAssertTrue(testCase, chooseTwo(64) == 2016);
    CuAssertTrue(testCase, chooseTwo(100) == 4950);
    CuAssertTrue(testCase, chooseTwo(1000) == 499500);
    CuAssertTrue(testCase, chooseTwo(4623824) == 10689871879576);
    uint64_t *cta = buildChooseTwoArray();
    for (uint32_t i = 0; i < 101; ++i) {
        CuAssertTrue(testCase, chooseTwo(i) == cta[i]);
    }
    // clean up
    free(cta);
}
static void checkPair(CuTest *testCase, stSortedSetIterator *sit, const char *a, const char *b, uint32_t i) {
    APair *p = NULL;
    p = stSortedSet_getNext(sit);
    CuAssertTrue(testCase, p != NULL);
    CuAssertTrue(testCase, strcmp(p->seq1, a) == 0);
    CuAssertTrue(testCase, strcmp(p->seq2, b) == 0);
    CuAssertTrue(testCase, p->pos1 == i);
    CuAssertTrue(testCase, p->pos2 == i);
}
static void createAndInsertPair(stSortedSet *pairs, const char *a, const char *b, uint32_t i) {
    APair *p = aPair_construct(stString_copy(a), stString_copy(b), i, i);
    stSortedSet_insert(pairs, aPair_copyConstruct(p));
    aPair_destruct(p);
}
static void test_pairSortComparison_0(CuTest *testCase) {
    stSortedSet *pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, (void(*)(void *)) aPair_destruct);
    APair *p = NULL;
    for (uint32_t i = 0; i < 5; ++i) {
        createAndInsertPair(pairs, "C", "D", i);
        createAndInsertPair(pairs, "B", "D", i);
        createAndInsertPair(pairs, "B", "C", i);
        createAndInsertPair(pairs, "A", "D", i);
        createAndInsertPair(pairs, "A", "C", i);
        createAndInsertPair(pairs, "A", "B", i);
        createAndInsertPair(pairs, ".", "z", i);
    }
    stSortedSetIterator *sit = stSortedSet_getIterator(pairs);
    for (uint32_t i = 0; i < 5; ++i) {
        checkPair(testCase, sit, ".", "z", i);
    }
    for (uint32_t i = 0; i < 5; ++i) {
        checkPair(testCase, sit, "A", "B", i);
        checkPair(testCase, sit, "A", "C", i);
        checkPair(testCase, sit, "A", "D", i);
    }
    for (uint32_t i = 0; i < 5; ++i) {
        checkPair(testCase, sit, "B", "C", i);
        checkPair(testCase, sit, "B", "D", i);
    }
    for (uint32_t i = 0; i < 5; ++i) {
        checkPair(testCase, sit, "C", "D", i);
    }
    stSortedSet_destructIterator(sit);
    APair *q = aPair_init();
    q->seq1 = stString_copy("A");
    q->pos1 = 1;
    p = stSortedSet_searchGreaterThanOrEqual(pairs, q);
    CuAssertTrue(testCase, p != NULL);
    sit = stSortedSet_getIteratorFrom(pairs, p);
    for (uint32_t i = 1; i < 5; ++i) {
        checkPair(testCase, sit, "A", "B", i);
        checkPair(testCase, sit, "A", "C", i);
        checkPair(testCase, sit, "A", "D", i);
    }
    aPair_destruct(q);
    stSortedSet_destructIterator(sit);
    q = aPair_init();
    q->seq1 = stString_copy("B");
    q->pos1 = 2;
    p = stSortedSet_searchGreaterThanOrEqual(pairs, q);
    CuAssertTrue(testCase, p != NULL);
    sit = stSortedSet_getIteratorFrom(pairs, p);
    for (uint32_t i = 2; i < 5; ++i) {
        checkPair(testCase, sit, "B", "C", i);
        checkPair(testCase, sit, "B", "D", i);
    }
    for (uint32_t i = 0; i < 5; ++i) {
        checkPair(testCase, sit, "C", "D", i);
    }
    aPair_destruct(q);
    stSortedSet_destructIterator(sit);
    q = aPair_init();
    q->seq1 = stString_copy("C");
    q->pos1 = 3;
    p = stSortedSet_searchGreaterThanOrEqual(pairs, q);
    CuAssertTrue(testCase, p != NULL);
    sit = stSortedSet_getIteratorFrom(pairs, p);
    for (uint32_t i = 3; i < 5; ++i) {
        checkPair(testCase, sit, "C", "D", i);
    }
    aPair_destruct(q);
    stSortedSet_destructIterator(sit);
    q = aPair_init();
    q->seq1 = stString_copy("D");
    q->pos1 = 3;
    p = stSortedSet_searchGreaterThanOrEqual(pairs, q);
    CuAssertTrue(testCase, p == NULL);
    aPair_destruct(q);
    // clean up
    
    stSortedSet_destruct(pairs);
    
}
CuSuite* comparatorAPI_TestSuite(void) {
    (void) (printMat);
    (void) (runningAverage);
    (void) (test_columnSampling_timing_0);
    (void) (test_mappingMatrixToArray_0);
    (void) (test_mappingArrayToMatrix_0);
    (void) (test_pairCounting_0);
    (void) (test_chooseTwoValues_0);
    (void) (test_columnSampling_timing_0);
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_mappingMatrixToArray_0);
    SUITE_ADD_TEST(suite, test_mappingArrayToMatrix_0);
    SUITE_ADD_TEST(suite, test_pairCounting_0);
    SUITE_ADD_TEST(suite, test_chooseTwoValues_0);
    SUITE_ADD_TEST(suite, test_pairSortComparison_0);
    // SUITE_ADD_TEST(suite, test_columnSampling_timing_0);
    return suite;
}
