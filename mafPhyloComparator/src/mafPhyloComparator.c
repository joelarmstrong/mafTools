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
#include "coalescences.h"
#include "blockTree.h"
#include "mafPhyloComparator.h"

void parseOpts(int argc, char **argv, PhyloOptions *opts);
void phyloOptions_destruct(PhyloOptions *opts);
void getLegitSequencesAndLengths(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash);

void parseOpts(int argc, char *argv[], PhyloOptions *opts) {
    st_logInfo("Parsing arguments\n");
    struct option longopts[] = {
        {"mafFile1", required_argument, NULL, 0},
        {"mafFile2", required_argument, NULL, 0},
        {"logLevel", required_argument, NULL, 0},
        {"numSamples", required_argument, NULL, 0},
        {"speciesTree", required_argument, NULL, 0},
        {"out", required_argument, NULL, 0},
        {"onlyLeaves", no_argument, NULL, 0},
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
            ret = sscanf(optarg, "%" PRIi64, &intArg);
            assert(ret == 1);
            opts->numSamples = intArg;
        } else if (strcmp(optName, "speciesTree") == 0) {
            opts->speciesTree = stTree_parseNewickString(optarg);
        } else if (strcmp(optName, "out") == 0) {
            opts->outFile = stString_copy(optarg);
        } else if (strcmp(optName, "onlyLeaves") == 0) {
            opts->onlyLeaves = true;
        }
    }
    if (opts->mafFile1 == NULL) {
        st_errAbort("--mafFile1 must be specified");
    }
    if (opts->mafFile2 == NULL) {
        st_errAbort("--mafFile2 must be specified");        
    }
    if (opts->speciesTree == NULL) {
        st_errAbort("--speciesTree is required (in newick format)");
    }
    if (opts->numSamples == 0) {
        opts->numSamples = 1000000;
    }
}

void phyloOptions_destruct(PhyloOptions *opts) {
    free(opts->mafFile1);
    free(opts->mafFile2);
    free(opts->outFile);
    stTree_destruct(opts->speciesTree);
    free(opts);
}

// Use the mafComparator API to build the legitSequences set and
// length hash -- but it expects a different options structure.
void getLegitSequencesAndLengths(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash) {
    Options *comparatorOpts = options_construct();
    comparatorOpts->mafFile1 = stString_copy(opts->mafFile1);
    comparatorOpts->mafFile2 = stString_copy(opts->mafFile2);
    buildSeqNamesSet(comparatorOpts, legitSequences, sequenceLengthHash);
    options_destruct(comparatorOpts);
}

// Check the species tree and exit and complain if it's invalid.
// All nodes must be labeled and there can be no duplicate names.
static void checkSpeciesTree(stTree *tree) {
    stSet *seen = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stList *nodes = stList_construct();
    fillListByReversePostOrder(tree, nodes, false);
    for (int64_t i = 0; i < stList_length(nodes); i++) {
        stTree *node = stList_get(nodes, i);
        if (stTree_getLabel(node) == NULL) {
            st_errAbort("All nodes in the species tree must be labeled.\n");
        }
        if (stSet_search(seen, (char *) stTree_getLabel(node))) {
            st_errAbort("Duplicate node found in species tree: %s", stTree_getLabel(node));
        }
        stSet_insert(seen, stString_copy(stTree_getLabel(node)));
    }

    stList_destruct(nodes);
    stSet_destruct(seen);
}

int main(int argc, char *argv[]) {
    PhyloOptions *opts = st_calloc(1, sizeof(PhyloOptions));
    parseOpts(argc, argv, opts);
    checkSpeciesTree(opts->speciesTree);

    // TODO: verify that the MAF has a tree for each block and the
    // tree is the format we need?
    st_logDebug("Getting legit sequences and lengths\n");
    stHash *sequenceLengthHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    stSet *legitSequences = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    getLegitSequencesAndLengths(opts, legitSequences, sequenceLengthHash);

    compareMAFCoalescences(opts, legitSequences, sequenceLengthHash, opts->onlyLeaves);

    // Clean up.
    phyloOptions_destruct(opts);
    stHash_destruct(sequenceLengthHash);
    stSet_destruct(legitSequences);
}
