#include <stdlib.h>
#include "CuTest.h"
#include "sonLib.h"
#include "sharedMaf.h"
#include "blockTree.h"
#include "test.blockTree.h"

static char *getRandomString(int64_t length) {
    char *ret = st_malloc(length + 1);
    for (int64_t i = 0; i < length; i++) {
        ret[i] = st_randomInt(0x41, 0x41 + 26);
    }
    ret[length] = '\0';
    return ret;
}

static stTree *getRandomTree(int64_t maxDepth) {
    stTree *tree = stTree_construct();
    stTree_setLabel(tree, getRandomString(8));
    if (maxDepth > 1) {
        int64_t numChildren = st_randomInt64(0, 10);
        for (int64_t i = 0; i < numChildren; i++) {
            stTree *child = getRandomTree(maxDepth - 1);
            stTree_setParent(child, tree);
        }
    }
    return tree;
}

static void test_parseTreeFromBlockStart_random(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 100; testNum++) {
        stTree *tree = getRandomTree(6);
        char *newick = stTree_getNewickTreeString(tree);
        char *headerLine = stString_print("a score=0.000 tree=\"%s\"", newick);
        mafLine_t *mafLine = maf_newMafLineFromString(headerLine, 0);
        
        // CuAssertStrEquals will run out of buffer space on large trees
        CuAssertTrue(testCase, strcmp(newick, parseTreeFromBlockStart(mafLine)) == 0);
        stTree_destruct(tree);
        free(newick);
        maf_destroyMafLineList(mafLine);
    }
}

static void test_parseTreeFromBlockStart(CuTest *testCase) {
    // check that tokens work correctly.
    mafLine_t *mafLine = maf_newMafLineFromString("a tree=\"(node1:0.1, node2:0.2);\"", 0);
    char *result = parseTreeFromBlockStart(mafLine);
    CuAssertStrEquals(testCase, "(node1:0.1, node2:0.2);", result);
    free(result);
    maf_destroyMafLineList(mafLine);

    // single token
    mafLine = maf_newMafLineFromString("a tree=\"(node1:0.1,node2:0.2);\"", 0);
    result = parseTreeFromBlockStart(mafLine);
    CuAssertStrEquals(testCase, "(node1:0.1,node2:0.2);", result);
    free(result);
    maf_destroyMafLineList(mafLine);
}

static void test_getMRCA(CuTest *testCase) {
    stTree *tree = stTree_parseNewickString("((a,b)c, (d, e)f, g)h;");
    stTree *node1 = stTree_findChild(tree, "a");
    stTree *node2 = stTree_findChild(tree, "b");
    CuAssertStrEquals(testCase, "c", stTree_getLabel(getMRCA(node1, node2)));

    node1 = stTree_findChild(tree, "b");
    node2 = stTree_findChild(tree, "a");
    CuAssertStrEquals(testCase, "c", stTree_getLabel(getMRCA(node1, node2)));

    node1 = stTree_findChild(tree, "a");
    node2 = stTree_findChild(tree, "c");
    CuAssertStrEquals(testCase, "c", stTree_getLabel(getMRCA(node1, node2)));

    node1 = stTree_findChild(tree, "a");
    node2 = stTree_findChild(tree, "a");
    CuAssertStrEquals(testCase, "a", stTree_getLabel(getMRCA(node1, node2)));

    node1 = stTree_findChild(tree, "a");
    node2 = stTree_findChild(tree, "d");
    CuAssertStrEquals(testCase, "h", stTree_getLabel(getMRCA(node1, node2)));

    node1 = stTree_findChild(tree, "a");
    node2 = stTree_findChild(tree, "g");
    CuAssertStrEquals(testCase, "h", stTree_getLabel(getMRCA(node1, node2)));

    stTree_destruct(tree);
}

static void test_getSeqToBlockRows_and_getNodeFromPosition(CuTest *testCase) {
    char *blockStr = "a tree=\"((a, a, b, c)d, (b, c)d)e;\"\n"
                     "s a 100 3 + 200 ggg\n"
                     "s a 197 3 - 200 ggg\n"
                     "s b 0 3 + 200 ggg\n"
                     "s c 0 3 + 200 ggg\n"
                     "s d 0 3 + 200 ggg\n"
                     "s b 100 3 + 200 ggg\n"
                     "s c 100 3 + 200 ggg\n"
                     "s d 100 3 + 200 ggg\n"
                     "s e 0 3 + 200 ggg\n";
    mafBlock_t *block = maf_newMafBlockFromString(blockStr, 0);
    stTree *tree = stTree_parseNewickString("((a,a,b,c)d, (b,c)d)e;");
    stHash *seqToBlockRows = getSeqToBlockRows(block, tree);

    // Row 1 -- a:100-103
    stTree *node = getNodeFromPosition(seqToBlockRows, "a", 100);
    CuAssertTrue(testCase, node == stTree_getChild(stTree_getChild(tree, 0), 0));
    stTree *prevNode = node;
    node = getNodeFromPosition(seqToBlockRows, "a", 101);
    CuAssertTrue(testCase, node == prevNode);
    node = getNodeFromPosition(seqToBlockRows, "a", 102);
    CuAssertTrue(testCase, node == prevNode);

    // Row 2 -- a:0-3
    node = getNodeFromPosition(seqToBlockRows, "a", 0);
    CuAssertTrue(testCase, node == stTree_getChild(stTree_getChild(tree, 0), 1));
    prevNode = node;
    node = getNodeFromPosition(seqToBlockRows, "a", 1);
    CuAssertTrue(testCase, node == prevNode);
    node = getNodeFromPosition(seqToBlockRows, "a", 2);
    CuAssertTrue(testCase, node == prevNode);

    // Row 3 -- b:0-3
    node = getNodeFromPosition(seqToBlockRows, "b", 0);
    CuAssertTrue(testCase, node == stTree_getChild(stTree_getChild(tree, 0), 2));
    prevNode = node;
    node = getNodeFromPosition(seqToBlockRows, "b", 1);
    CuAssertTrue(testCase, node == prevNode);
    node = getNodeFromPosition(seqToBlockRows, "b", 2);
    CuAssertTrue(testCase, node == prevNode);

    // Row 7 -- c:100-103
    node = getNodeFromPosition(seqToBlockRows, "c", 100);
    CuAssertTrue(testCase, node == stTree_getChild(stTree_getChild(tree, 1), 1));
    prevNode = node;
    node = getNodeFromPosition(seqToBlockRows, "c", 101);
    CuAssertTrue(testCase, node == prevNode);
    node = getNodeFromPosition(seqToBlockRows, "c", 102);
    CuAssertTrue(testCase, node == prevNode);

    // Row 8 -- d:100-103
    node = getNodeFromPosition(seqToBlockRows, "d", 100);
    CuAssertTrue(testCase, node == stTree_getChild(tree, 1));
    prevNode = node;
    node = getNodeFromPosition(seqToBlockRows, "d", 101);
    CuAssertTrue(testCase, node == prevNode);
    node = getNodeFromPosition(seqToBlockRows, "d", 102);
    CuAssertTrue(testCase, node == prevNode);

    // Row 9 -- e:0-3
    node = getNodeFromPosition(seqToBlockRows, "e", 0);
    CuAssertTrue(testCase, node == tree);
    prevNode = node;
    node = getNodeFromPosition(seqToBlockRows, "e", 1);
    CuAssertTrue(testCase, node == prevNode);
    node = getNodeFromPosition(seqToBlockRows, "e", 2);
    CuAssertTrue(testCase, node == prevNode);

    // Not in block = NULL
    node = getNodeFromPosition(seqToBlockRows, "e", 50);
    CuAssertTrue(testCase, node == NULL);

    maf_destroyMafBlockList(block);
    stTree_destruct(tree);
    stHash_destruct(seqToBlockRows);
}

static void test_buildNameToNodeHash_random(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 100; testNum++) {
        stTree *tree = getRandomTree(3);
        stHash *nameToNode = buildNameToNodeHash(tree);
        stHashIterator *hashIt = stHash_getIterator(nameToNode);
        char *name;
        while ((name = stHash_getNext(hashIt)) != NULL) {
            stTree *node = stHash_search(nameToNode, name);
            CuAssertStrEquals(testCase, name, stTree_getLabel(node));
        }

        stHash_destruct(nameToNode);
        stTree_destruct(tree);
    }
}

static void test_isAncestor(CuTest *testCase) {
    stTree *tree = stTree_parseNewickString("(a,b,(c,d)e)f;");
    stHash *nameToNode = buildNameToNodeHash(tree);
    CuAssertTrue(testCase, isAncestor("f", "a", nameToNode));
    CuAssertTrue(testCase, isAncestor("f", "d", nameToNode));
    CuAssertTrue(testCase, isAncestor("e", "d", nameToNode));
    CuAssertTrue(testCase, !isAncestor("a", "d", nameToNode));
    CuAssertTrue(testCase, !isAncestor("a", "f", nameToNode));
    stHash_destruct(nameToNode);
    stTree_destruct(tree);
}

CuSuite *blockTree_TestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_parseTreeFromBlockStart);
    SUITE_ADD_TEST(suite, test_parseTreeFromBlockStart_random);
    SUITE_ADD_TEST(suite, test_getMRCA);
    SUITE_ADD_TEST(suite, test_getSeqToBlockRows_and_getNodeFromPosition);
    SUITE_ADD_TEST(suite, test_buildNameToNodeHash_random);
    SUITE_ADD_TEST(suite, test_isAncestor);
    return suite;
}
