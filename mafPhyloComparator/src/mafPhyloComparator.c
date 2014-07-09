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

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include "sonLib.h"
#include "common.h"
#include "sharedMaf.h"
#include "comparatorAPI.h"
#include "mafPhyloComparator.h"

void parseOpts(int argc, char **argv, PhyloOptions *opts);
char *parseTreeFromBlockStart(mafLine_t *line);
void sampleCoalescences(char *mafFileName, stSet *coalescences, double acceptProbability, stSet *legitSequences, stHash *sequenceLengthHash);
void getLegitSequencesAndLengths(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash);

void parseOpts(int argc, char *argv[], PhyloOptions *opts) {
    st_logInfo("Parsing arguments\n");
    struct option longopts[] = {
        {"mafFile1", required_argument, NULL, 0},
        {"mafFile2", required_argument, NULL, 0},
        {"logLevel", required_argument, NULL, 0},
        {"numSamples", required_argument, NULL, 0},
        {0, 0, 0, 0}
    };
    int longindex;
    while (getopt_long(argc, argv, "", longopts, &longindex) != -1) {
        const char *optName = longopts[longindex].name;
        if (strcmp(optName, "mafFile1") == 0) {
            opts->mafFile1 = stString_copy(optarg);
        } else if (strcmp(optName, "mafFile2") == 0) {
            opts->mafFile2 = stString_copy(optarg);
        } else if (strcmp(optName, "logLevel") == 0) {
            st_setLogLevelFromString(optarg);
        } else if (strcmp(optName, "numSamples") == 0) {
            int64_t intArg;
            int ret;
            ret = sscanf("%" PRIi64, optarg, &intArg);
            assert(ret == 1);
            opts->numSamples = intArg;
        }
    }
    if (opts->mafFile1 == NULL) {
        st_errAbort("--mafFile1 must be specified");
    }
    if (opts->mafFile2 == NULL) {
        st_errAbort("--mafFile2 must be specified");        
    }
    if (opts->numSamples == 0) {
        opts->numSamples = 1000000;
    }
}

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
            stList_append(treeTokens, stString_copy(tok + 6));
        } else if (inTreeString) {
            int64_t len = strlen(tok);
            if (tok[len - 1] == '"') {
                inTreeString = false;
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
    
}

// Get a set of coalescences from a set of APairs given a gene tree.
void coalescencesFromPairs(stTree *tree, stSortedSet *pairs, stHash *intervalToNode, stSet *coalescences) {

    stSortedSet *pairsIt = stSortedSet_getIterator(pairs);    
}

// For use in a set containing rows from only one sequence, so
// comparing sequence strings is not necessary.
int BlockRow_cmp(BlockRow *row1, BlockRow *row2) {
    assert(strcmp(row1->seq, row2->seq) == 0);
    if (row2->start > row1->start) {
        assert(row2->start >= row1->end); // Intervals should never overlap.
        assert(row2->end > row1->end);
        return 1;
    } else if (row2->start == row1->start) {
        return 0;
    } else {
        assert(row1->start >= row2->end); // Intervals should never overlap.
        assert(row1->end > row2->end);
        return -1;
    }
}

stTree *getNodeFromPosition(stHash *seqToBlockRows, const char *seq, uint64_t pos) {
    stSortedSet *blockRows = stHash_search(seqToBlockRows, seq);
    if (blockRows == NULL) {
        return NULL;
    }

    // Need to allocate a fake interval to search for an interval
    // containing the position. Pretty hacky.
    BlockRow *tmp = st_malloc(1, sizeof(BlockRow));
    tmp->seq = seq;
    tmp->start = pos;
    tmp->end = pos + 1;
    BlockRow *blockRow = stSortedSet_searchLessThanOrEqual(blockRows, tmp);
    free(tmp);
    if (blockRow->end > pos) {
        // Found a block row with an interval containing this position.
        return blockRow;
    } else {
        return NULL;
    }
}

// Get a hash from seq -> stSortedSet that contains BlockRow
// structs. (kind of a poor man's map so that intervals can be mapped
// to tree nodes.)
stHash *getBlockRows(mafBlock_t *block, stTree *tree) {
    stHash *ret = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void *)) stSortedSet_destruct);
    // Walk through the block and assign rows to nodes in a post-order fashion.
    
    return ret;
}

// Walk through the given maf file, sampling pairs and recording where
// in the tree they coalesce.
void sampleCoalescences(char *mafFileName, stSet *coalescences, double acceptProbability, stSet *legitSequences, stHash *sequenceLengthHash) {
    mafFileApi_t *mafFile = maf_newMfa(mafFileName, "r");
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    mafBlock_t *block;
    while ((block = maf_readBlock(mafFile)) != NULL) {
        // Parse out tree header
        mafLine_t *line = maf_mafBlock_getHeadLine(block);
        if(maf_mafLine_getType(line) == 'h') {
            // header line, we should just skip this block
            maf_destroyMafBlockList(block);
            continue;
        }
        assert(maf_mafLine_getType(line) == 'a');

        // Tree stuff
        char *newickString = parseTreeFromBlockStart(line);
        st_logDebug("Got gene tree %s from block\n", newickString);
        stTree *tree = stTree_parseNewickString(newickString);
        stHash *intervalToNode = getIntervalToNode(block, tree);

        // Use existing mafComparator api to get pairs
        stSortedSet *pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, (void(*)(void *)) aPair_destruct);
        uint64_t numPairs = 0;
        walkBlockSamplingPairs(mafFileName, block, pairs, acceptProbability, legitSequences, chooseTwoArray, &numPairs, sequenceLengthHash);
        st_logDebug("Sampled %d pairs from block\n", numPairs);

        coalescencesFromPairs(tree, pairs, intervalToNode, coalescences);

        stSortedSet_destruct(pairs);
        maf_destroyMafBlockList(block);
    }
    maf_destroyMfa(mafFile);
}

// Use the mafComparator API to build the legitSequences set and
// length hash -- but it expects a different options structure.
void getLegitSequencesAndLengths(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash) {
    Options *comparatorOpts = options_construct();
    comparatorOpts->mafFile1 = opts->mafFile1;
    comparatorOpts->mafFile2 = opts->mafFile2;
    buildSeqNamesSet(comparatorOpts, legitSequences, sequenceLengthHash);
    options_destruct(comparatorOpts);
}

int main(int argc, char *argv[]) {
    PhyloOptions *opts = st_calloc(1, sizeof(PhyloOptions));
    parseOpts(argc, argv, opts);

    // TODO: verify that the MAF has a tree for each block and the
    // tree is the format we need?
    stHash *sequenceLengthHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    stSet *legitSequences = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    getLegitSequencesAndLengths(opts, legitSequences, sequenceLengthHash);

    // Sample coalescences from the MAF (a pair of sequences from a
    // block and what genome they coalesce in)
    st_logInfo("Sampling coalescences\n");
    stSet *coalescences = stSet_construct();
    double acceptProbability = opts->numSamples / countPairsInMaf(opts->mafFile1, legitSequences);
    sampleCoalescences(opts->mafFile1, coalescences, acceptProbability, legitSequences, sequenceLengthHash);

    /* findMatchingCoalescences(mafFile2, coalescences, matchedCoalescences); */

    /* scoreMatchedCoalescences(coalescences, matchedCoalescences); */

    free(opts);
}
