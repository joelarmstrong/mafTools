/*
 * Copyright (C) 2009-2014 by
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
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
#include <stdlib.h>
#include "sonLib.h"
#include "sharedMaf.h"
#include "common.h"
#include "blockTree.h"

// Get the gene tree from the block start. If there is no species
// tree in the line, raise an error.
// Trees (at least for now) must be quoted.
char *parseTreeFromBlockStart(mafLine_t *line) {
    // FIXME: super sloppy
    char *lineString = maf_mafLine_getLine(line);
    char **linePtr = &lineString;
    char *tok = de_strtok(linePtr, ' ');
    bool inTreeString = false, seenTree = false;
    stList *treeTokens = stList_construct3(0, free);
    while (tok != NULL) {
        if (strncmp(tok, "tree=", 5) == 0) {
            // Found a tree start
            seenTree = true;
            if (tok[5] != '"') {
                st_errAbort("Maf block at line %d contains a tree parameter, but it is not enclosed in quotes.\n", maf_mafLine_getLineNumber(line));
            }
            inTreeString = true;

            int64_t len = strlen(tok);
            if (tok[len - 1] == '"') {
                inTreeString = false;
                tok [len - 1] = '\0';
            }
            stList_append(treeTokens, stString_copy(tok + 6));
        } else if (inTreeString) {
            int64_t len = strlen(tok);
            if (tok[len - 1] == '"') {
                inTreeString = false;
                tok[len - 1] = '\0';
            }
            stList_append(treeTokens, stString_copy(tok));
        }
        free(tok);
        tok = de_strtok(linePtr, ' ');
    }
    if (inTreeString) {
        st_errAbort("Maf block start at line %d has a partially "
                    "quoted tree string.\n",
                    maf_mafLine_getLineNumber(line));
    }
    if (!seenTree) {
        st_errAbort("Maf block start at line %d does not contain a "
                    "tree. mafPhyloComparator requires a tree for "
                    "all blocks.\n", maf_mafLine_getLineNumber(line));
    }

    char *newick = stString_join2(" ", treeTokens);
    stList_destruct(treeTokens);
    return newick;
}

// Will be super slow. TODO: speed up by including useful info in the
// clientdata field.
stTree *getMRCA(stTree *node1, stTree *node2) {
    stTree *ret = NULL;

    // Get all of node1's parents
    stSet *parents = stSet_construct();
    stTree *node = node1;
    do {
        stSet_insert(parents, node);
    } while ((node = stTree_getParent(node)) != NULL);

    // Compare node2's parents against node1's parents and find the
    // first match
    node = node2;
    do {
        if (stSet_search(parents, node) != NULL) {
            ret = node;
            break;
        }
    } while ((node = stTree_getParent(node)) != NULL);

    stSet_destruct(parents);
    return ret;
}

// For use in a set containing rows from only one sequence, so
// comparing sequence strings is not necessary.
int BlockRow_cmp(const BlockRow *row1, const BlockRow *row2) {
    assert(strcmp(row1->species, row2->species) == 0);
    if (row2->start > row1->start) {
        return -1;
    } else if (row2->start == row1->start) {
        return 0;
    } else {
        return 1;
    }
}

void BlockRow_destruct(BlockRow *row) {
    free(row->species);
    free(row);
}

stTree *getNodeFromPosition(stHash *seqToBlockRows, const char *seq, uint64_t pos) {
    stSortedSet *blockRows = stHash_search(seqToBlockRows, (void *) seq);
    if (blockRows == NULL) {
        return NULL;
    }

    // Need to allocate a fake interval to search for an interval
    // containing the position. Pretty hacky.
    BlockRow *tmp = st_malloc(sizeof(BlockRow));
    tmp->species = stString_copy(seq);
    tmp->start = pos;
    tmp->end = pos + 1;
    BlockRow *blockRow = stSortedSet_searchLessThanOrEqual(blockRows, tmp);
    free(tmp->species);
    free(tmp);
    if (blockRow == NULL) {
        return NULL;
    }
    if (blockRow->end > pos) {
        // Found a block row with an interval containing this position.
        return blockRow->node;
    } else {
        // The closest interval doesn't contain this position.
        return NULL;
    }
}

// fill by reversed post-order (i.e. pre-order with the order we visit
// the children in reversed), so we can pop off the end of the list
// rather than the start.
void fillListByReversePostOrder(stTree *tree, stList *list) {
    stList_append(list, tree);
    for (int64_t i = stTree_getChildNumber(tree) - 1; i >= 0; i--) {
        fillListByReversePostOrder(stTree_getChild(tree, i), list);
    }
}

// Get a hash from seq -> stSortedSet that contains BlockRow
// structs. (kind of a poor man's map so that intervals can be mapped
// to tree nodes.)
// The tree should correspond to the rows in the block in post-order.
stHash *getSeqToBlockRows(mafBlock_t *block, stTree *tree) {
    stHash *ret = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) stSortedSet_destruct);
    // Walk through the block and assign rows to nodes in a post-order fashion.
    stList *nodes = stList_construct();
    // Need to pop off the end of the list for efficiency reasons.
    fillListByReversePostOrder(tree, nodes);
    st_logDebug("Found %" PRIi64 " species in the block tree\n", stList_length(nodes));
    mafLine_t *line = maf_mafBlock_getHeadLine(block);
    for (;;) {
        st_logDebug("Processing line %" PRIi64 "\n", maf_mafLine_getLineNumber(line));
        if (maf_mafLine_getType(line) != 's') {
            // Not a sequence line; skip.
            line = maf_mafLine_getNext(line);
            continue;
        }

        stTree *node = stList_pop(nodes);
        char *header = maf_mafLine_getSpecies(line);
        stSortedSet *sortedSet = stHash_search(ret, header);
        if (sortedSet == NULL) {
            // Initalize hash entry for this header.
            sortedSet = stSortedSet_construct3((int (*)(const void *, const void *)) BlockRow_cmp, (void (*)(void *)) BlockRow_destruct);
            stHash_insert(ret, stString_copy(header), sortedSet);
        }

        BlockRow *blockRow = st_malloc(sizeof(BlockRow));
        blockRow->species = stString_copy(header);
        if (maf_mafLine_getStrand(line) == '+') {
            blockRow->start = maf_mafLine_getStart(line);
            blockRow->end = blockRow->start + maf_mafLine_getLength(line);
        } else {
            // Reversed
            blockRow->start = maf_mafLine_getSourceLength(line) - maf_mafLine_getStart(line) - maf_mafLine_getLength(line);
            blockRow->end = blockRow->start + maf_mafLine_getLength(line);
        }
        assert(blockRow->start >= 0);
        assert(blockRow->start < blockRow->end);
        assert(strcmp(blockRow->species, stTree_getLabel(node)) == 0);
        blockRow->node = node;
        stSortedSet_insert(sortedSet, blockRow);

        if (line == maf_mafBlock_getTailLine(block)) {
            break;
        } else {
            line = maf_mafLine_getNext(line);
        }
    }
    assert(stList_length(nodes) == 0); // The # of nodes should = the
                                       // # of rows
    stList_destruct(nodes);
    return ret;
}

// Get a mapping from tree label to node.
stHash *buildNameToNodeHash(stTree *tree) {
    stHash *nameToNode = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    stList *nodes = stList_construct();
    fillListByReversePostOrder(tree, nodes);
    for (int64_t i = 0; i < stList_length(nodes); i++) {
        stTree *node = stList_get(nodes, i);
        assert(stHash_search(nameToNode, (char *) stTree_getLabel(node)) == NULL);
        stHash_insert(nameToNode, stString_copy(stTree_getLabel(node)), node);
    }
    stList_destruct(nodes);
    return nameToNode;
}

// returns true if name1 is an ancestor of name2. nameToNode is a hash
// from strings to their corresponding nodes in a species tree.
bool isAncestor(char *name1, char *name2, stHash *nameToNode) {
    stTree *tree1 = stHash_search(nameToNode, name1);
    assert(tree1 != NULL);
    stTree *tree2 = stHash_search(nameToNode, name2);
    assert(tree2 != NULL);

    if (tree1 == tree2) {
        return false;
    } else {
        return getMRCA(tree1, tree2) == tree1;
    }
}
