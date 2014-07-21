#include "CuTest.h"
#include "sonLib.h"
#include "sharedMaf.h"
#include "comparatorAPI.h"
#include "blockTree.h"
#include "coalescences.h"
#include "test.coalescences.h"

static void test_coalescencesFromPairs(CuTest *testCase) {
    // get the basic structures (block tree, seqToBlockRow hash, etc)
    const char *blockStr = "a tree=\"((a, b)d, (b, c)d)e;\"\n"
                           "s a 100 3 + 200 ggg\n"
                           "s b 0 3 + 200 ggg\n"
                           "s d 0 3 + 200 ggg\n"
                           "s b 100 3 + 200 ggg\n"
                           "s c 100 3 + 200 ggg\n"
                           "s d 100 3 + 200 ggg\n"
                           "s e 0 3 + 200 ggg\n";
    mafBlock_t *block = maf_newMafBlockFromString(blockStr, 0);
    stTree *tree = stTree_parseNewickString("((a, b)d,(b, c)d)e;");
    stHash *seqToBlockRows = getSeqToBlockRows(block, tree, false);

    // Construct the test pairs
    stSortedSet *pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, (void(*)(void *)) aPair_destruct);
    APair *pair = aPair_init();
    aPair_fillOut(pair, stString_copy("a"), stString_copy("b"), 100, 1);
    stSortedSet_insert(pairs, pair);
    pair = aPair_init();
    aPair_fillOut(pair, stString_copy("a"), stString_copy("b"), 100, 100);
    stSortedSet_insert(pairs, pair);

    // Get coalescences and check that they are correct
    stSortedSet *coalescences = stSortedSet_construct3((int (*)(const void *, const void *)) coalescence_cmp, (void (*)(void *)) coalescence_destruct);
    coalescencesFromPairs(pairs, seqToBlockRows, coalescences);
    CuAssertIntEquals(testCase, 2, stSortedSet_size(coalescences));
    // First coalescence -- a:100,b:1. Should coalesce at species "d"
    Coalescence *coal = stSortedSet_getFirst(coalescences);
    CuAssertStrEquals(testCase, "a", coal->seq1);
    CuAssertIntEquals(testCase, 100, coal->pos1);
    CuAssertStrEquals(testCase, "b", coal->seq2);
    CuAssertIntEquals(testCase, 1, coal->pos2);
    CuAssertStrEquals(testCase, "d", coal->mrca);
    // Second coalescence -- a:100,b:100. Should coalesce at species "e"
    coal = stSortedSet_getLast(coalescences);
    CuAssertStrEquals(testCase, "a", coal->seq1);
    CuAssertIntEquals(testCase, 100, coal->pos1);
    CuAssertStrEquals(testCase, "b", coal->seq2);
    CuAssertIntEquals(testCase, 100, coal->pos2);
    CuAssertStrEquals(testCase, "e", coal->mrca);

    // Clean up
    stSortedSet_destruct(coalescences);
    stSortedSet_destruct(pairs);
    stTree_destruct(tree);
    maf_destroyMafBlockList(block);
}

CuSuite *coalescences_TestSuite(void) {
    CuSuite *ret = CuSuiteNew();
    SUITE_ADD_TEST(ret, test_coalescencesFromPairs);

    return ret;
}
