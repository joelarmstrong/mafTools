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
#include "common.h"
#include "sharedMaf.h"
#include "comparatorAPI.h"
#include "blockTree.h"
#include "mafPhyloComparator.h"
#include "coalescences.h"

static void sampleCoalescences(char *mafFileName, stSortedSet *coalescences, double acceptProbability, stSet *legitSequences, stHash *sequenceLengthHash, bool onlyLeaves);
static stSortedSet *sortedSetFromSet(stSet *set);
static stSortedSet *pairsFromCoalescences(stSortedSet *coalescences);
static stSortedSet *findMatchingCoalescences(char *mafFileName, stSortedSet *coalescences, stSet *legitSequences, bool onlyLeaves);
static CoalResult *coalResult_init(const char *seq);
static void coalResult_destruct(CoalResult *coalResult);
static CoalResult *getResultForSequence(stHash *seqResultsHash, char *seq);
static void buildCoalescenceResults(stSortedSet *sampledCoalescences, stSortedSet *matchedCoalescences, stTree *speciesTree, CoalResult *aggregateResults, stHash *seqResults);
static void reportCoalescenceResult(const char *tag, CoalResult *result, FILE *f);
static void reportCoalescenceResults(CoalResult *aggregateResults, stHash *seqResults, const char *mafFile1, const char *mafFile2, FILE *f);

// Get a set of coalescences from a sorted set of APairs given a gene
// tree.
void coalescencesFromPairs(stSortedSet *pairs, stHash *seqToBlockRows, stSortedSet *coalescences) {
    stSortedSetIterator *pairsIt = stSortedSet_getIterator(pairs);
    APair *pair;
    while ((pair = stSortedSet_getNext(pairsIt)) != NULL) {
        Coalescence *coalescence = st_malloc(sizeof(Coalescence));
        stTree *node1 = getNodeFromPosition(seqToBlockRows, pair->seq1, pair->pos1);
        stTree *node2 = getNodeFromPosition(seqToBlockRows, pair->seq2, pair->pos2);
        stTree *mrca = getMRCA(node1, node2);
        coalescence->seq1 = stString_copy(pair->seq1);
        coalescence->pos1 = pair->pos1;
        coalescence->seq2 = stString_copy(pair->seq2);
        coalescence->pos2 = pair->pos2;
        coalescence->mrca = stString_copy(stTree_getLabel(mrca));
        stSortedSet_insert(coalescences, coalescence);
    }

    stSortedSet_destructIterator(pairsIt);
}

// Sample coalescences from a block.
void walkBlockSamplingCoalescences(char *mafFileName, mafBlock_t *block, stSortedSet *coalescences, double acceptProbability, stSet *legitSequences, stHash *sequenceLengthHash, uint64_t *chooseTwoArray, bool onlyLeaves) {
    // Parse out tree header
    mafLine_t *line = maf_mafBlock_getHeadLine(block);
    assert(maf_mafLine_getType(line) == 'a');
    char *newickString = parseTreeFromBlockStart(line);
    st_logDebug("Got gene tree %s from block\n", newickString);
    stTree *tree = stTree_parseNewickString(newickString);
    stHash *seqToBlockRows = getSeqToBlockRows(block, tree, onlyLeaves);

    // Use existing mafComparator api to get pairs
    stSortedSet *pairs = stSortedSet_construct3((int(*)(const void *, const void *)) aPair_cmpFunction, (void(*)(void *)) aPair_destruct);
    uint64_t numPairs = 0;
    walkBlockSamplingPairs(mafFileName, block, pairs, acceptProbability, legitSequences, chooseTwoArray, &numPairs, sequenceLengthHash);
    st_logDebug("Sampled %" PRIi64 " of %" PRIi64 " pairs from block\n", stSortedSet_size(pairs), numPairs);

    coalescencesFromPairs(pairs, seqToBlockRows, coalescences);

    free(newickString);
    stTree_destruct(tree);
    stHash_destruct(seqToBlockRows);
    stSortedSet_destruct(pairs);    
}

// Walk through the given maf file, sampling pairs and recording where
// in the tree they coalesce.
static void sampleCoalescences(char *mafFileName, stSortedSet *coalescences, double acceptProbability, stSet *legitSequences, stHash *sequenceLengthHash, bool onlyLeaves) {
    mafFileApi_t *mafFile = maf_newMfa(mafFileName, "r");
    uint64_t *chooseTwoArray = buildChooseTwoArray();
    mafBlock_t *block;
    while ((block = maf_readBlock(mafFile)) != NULL) {
        mafLine_t *line = maf_mafBlock_getHeadLine(block);
        if (maf_mafLine_getType(line) == 'h') {
            // header line, we should just skip this block
            maf_destroyMafBlockList(block);
            continue;
        }

        walkBlockSamplingCoalescences(mafFileName, block, coalescences, acceptProbability, legitSequences, sequenceLengthHash, chooseTwoArray, onlyLeaves);
        maf_destroyMafBlockList(block);
    }

    free(chooseTwoArray);
    maf_destroyMfa(mafFile);
}

// Convenience function to get a sorted set (keyed by pointer) from a
// stSet. Not ideal to use very often since it is O(n).
static stSortedSet *sortedSetFromSet(stSet *set) {
    stSetIterator *setIt = stSet_getIterator(set);
    stSortedSet *ret = stSortedSet_construct();

    void *item;
    while((item = stSet_getNext(setIt)) != NULL) {
        stSortedSet_insert(ret, item);
    }

    stSet_destructIterator(setIt);
    return ret;
}

static stSortedSet *pairsFromCoalescences(stSortedSet *coalescences) {
    stSortedSet *ret = stSortedSet_construct3((int (*)(const void *, const void *))aPair_cmpFunction, (void (*)(void *)) aPair_destruct);
    stSortedSetIterator *setIt = stSortedSet_getIterator(coalescences);
    Coalescence *coal;
    while((coal = stSortedSet_getNext(setIt)) != NULL) {
        APair *pair = aPair_init();
        aPair_fillOut(pair, stString_copy(coal->seq1), stString_copy(coal->seq2), coal->pos1, coal->pos2);
        stSortedSet_insert(ret, pair);
    }

    stSortedSet_destructIterator(setIt);
    return ret;
}

static stSortedSet *findMatchingCoalescences(char *mafFileName, stSortedSet *coalescences, stSet *legitSequences, bool onlyLeaves) {
    stSortedSet *matchingCoalescences = stSortedSet_construct3((int (*)(const void *, const void *)) coalescence_cmp, (void (*)(void *)) coalescence_destruct);

    // Get aPairs from the sampled coalescences.
    stSortedSet *pairs = pairsFromCoalescences(coalescences);
    st_logDebug("Converted %" PRIi64 " coalescences back to pairs\n", stSortedSet_size(pairs));

    mafFileApi_t *mafFile = maf_newMfa(mafFileName, "r");
    mafBlock_t *block;
    while ((block = maf_readBlock(mafFile)) != NULL) {
        mafLine_t *line = maf_mafBlock_getHeadLine(block);
        if (maf_mafLine_getType(line) != 'a') {
            // Only looking for alignment blocks.
            maf_destroyMafBlockList(block);
            continue;
        }

        // Use existing mafComparator API to get matching pairs.
        stSet *matchingBlockPairs = stSet_construct();
        walkBlockTestingHomology(block, pairs, matchingBlockPairs,
                                 legitSequences, 0);
        st_logDebug("Got %" PRIi64 " matching pairs from the block\n", stSet_size(matchingBlockPairs));

        // Get the tree
        char *newickString = parseTreeFromBlockStart(line);
        stTree *blockTree = stTree_parseNewickString(newickString);
        stHash *seqToBlockRows = getSeqToBlockRows(block, blockTree, onlyLeaves);

        // Change from a set of positive pairs to a sorted set. NB:
        // this isn't really necessary, it's just to make the types
        // match, so if this is appreciably slow just duplicate the
        // coalescencesFromPairs into sorted set and unsorted set
        // versions.
        stSortedSet *matchingBlockPairsSorted = sortedSetFromSet(matchingBlockPairs);

        coalescencesFromPairs(matchingBlockPairsSorted,
                              seqToBlockRows, matchingCoalescences);

        free(newickString);
        stTree_destruct(blockTree);
        stHash_destruct(seqToBlockRows);
        stSet_destruct(matchingBlockPairs);
        stSortedSet_destruct(matchingBlockPairsSorted);
        maf_destroyMafBlockList(block);
    }

    stSortedSet_destruct(pairs);
    maf_destroyMfa(mafFile);

    return matchingCoalescences;
}

// Compare two coalescences, by comparing their first pairs, then
// their second pairs, then their MRCAs.
int coalescence_cmp(const Coalescence *coal1, const Coalescence *coal2) {
    int ret = strcmp(coal1->seq1, coal2->seq1);
    if (ret) {
        return ret;
    }
    if (coal1->pos1 > coal2->pos1) {
        return 1;
    } else if (coal1->pos1 < coal2->pos1) {
        return -1;
    }

    ret = strcmp(coal1->seq2, coal2->seq2);
    if (ret) {
        return ret;
    }
    if (coal1->pos2 > coal2->pos2) {
        return 1;
    } else if (coal1->pos2 < coal2->pos2) {
        return -1;
    }

    return strcmp(coal1->mrca, coal2->mrca);
}

void coalescence_destruct(Coalescence *coal) {
    free(coal->seq1);
    free(coal->seq2);
    free(coal->mrca);
    free(coal);
}

static CoalResult *coalResult_init(const char *seq) {
    CoalResult *ret = st_calloc(1, sizeof(CoalResult));
    ret->seq = stString_copy(seq);
    return ret;
}

static void coalResult_destruct(CoalResult *coalResult) {
    free(coalResult->seq);
    free(coalResult);
}

static CoalResult *getResultForSequence(stHash *seqResultsHash, char *seq) {
    CoalResult *ret = stHash_search(seqResultsHash, seq);
    if (ret == NULL) {
        // initialize results structure for this sequence.
        ret = coalResult_init(seq);
        stHash_insert(seqResultsHash, stString_copy(seq), ret);
    }
    return ret;
}

// Fill in result structures on a per-sequence and overall basis.
static void buildCoalescenceResults(stSortedSet *sampledCoalescences, stSortedSet *matchedCoalescences, stTree *speciesTree, CoalResult *aggregateResults, stHash *seqResults) {
    stHash *nameToSpecies = buildNameToNodeHash(speciesTree);

    // Go through each sampled coalescence (from maf A) and find the
    // coalescence from maf B that has the same seqs & positions, if
    // any. Then compare the genomes that they coalesce in to see if
    // they match.
    stSortedSetIterator *coalIt = stSortedSet_getIterator(sampledCoalescences);
    Coalescence *sampledCoal;
    while ((sampledCoal = stSortedSet_getNext(coalIt)) != NULL) {
        CoalResult *seq1Result = getResultForSequence(seqResults, sampledCoal->seq1);
        CoalResult *seq2Result = getResultForSequence(seqResults, sampledCoal->seq2);
        assert(strcmp(seq1Result->seq, sampledCoal->seq1) == 0);
        assert(strcmp(seq2Result->seq, sampledCoal->seq2) == 0);
        aggregateResults->sampledPairs++;
        seq1Result->sampledPairs++;
        seq2Result->sampledPairs++;
        Coalescence *matchedCoal = stSortedSet_search(matchedCoalescences, sampledCoal);
        if (matchedCoal != NULL) {
            // There is a matching pair in other MAF. Record
            // where the pairs coalesce in the two mafs.
            assert(strcmp(matchedCoal->seq1, sampledCoal->seq1) == 0);
            assert(strcmp(matchedCoal->seq2, sampledCoal->seq2) == 0);
            seq1Result->matchingCoalescences++;
            seq2Result->matchingCoalescences++;
            aggregateResults->matchingCoalescences++;
            if (strcmp(sampledCoal->mrca, matchedCoal->mrca) == 0) {
                seq1Result->identicalCoalescences++;
                seq2Result->identicalCoalescences++;
                aggregateResults->identicalCoalescences++;
            } else if (isAncestor(sampledCoal->mrca, matchedCoal->mrca, nameToSpecies)) {
                // Maf B's coalescence is earlier than maf A's.
                seq1Result->earlyCoalescences++;
                seq2Result->earlyCoalescences++;
                aggregateResults->earlyCoalescences++;
            } else {
                // Maf B's coalescence is later than maf A's.
                assert(isAncestor(matchedCoal->mrca, sampledCoal->mrca, nameToSpecies));
                seq1Result->lateCoalescences++;
                seq2Result->lateCoalescences++;
                aggregateResults->lateCoalescences++;
            }
        }
    }

    stSortedSet_destructIterator(coalIt);
    stHash_destruct(nameToSpecies);
}

static void reportCoalescenceResult(const char *tag, CoalResult *result, FILE *f) {
    fprintf(f, "<%s seq=\"%s\" sampled=\"%" PRIu64 "\" matching=\"%" PRIu64 "\" identicalCoalescences=\"%" PRIu64 "\" earlyCoalescences=\"%" PRIu64 "\" lateCoalescences=\"%" PRIu64 "\" avgIdentical=\"%f\" avgEarly=\"%f\" avgLate=\"%f\">\n", tag, result->seq, result->sampledPairs, result->matchingCoalescences, result->identicalCoalescences, result->earlyCoalescences, result->lateCoalescences, ((float)result->identicalCoalescences)/result->matchingCoalescences, ((float)result->earlyCoalescences)/result->matchingCoalescences, ((float)result->lateCoalescences)/result->matchingCoalescences);
}

static void reportCoalescenceResults(CoalResult *aggregateResults, stHash *seqResults, const char *mafFile1, const char *mafFile2, FILE *f) {
    assert(strcmp(aggregateResults->seq, "aggregate") == 0);
    fprintf(f, "<coalescenceTest mafFile1=\"%s\" mafFile2=\"%s\">\n", mafFile1, mafFile2);
    reportCoalescenceResult("aggregateCoalescenceResults", aggregateResults, f);
    stHashIterator *seqIt = stHash_getIterator(seqResults);
    char *seq;
    while ((seq = stHash_getNext(seqIt)) != NULL) {
        CoalResult *seqResult = stHash_search(seqResults, seq);
        assert(seqResult != NULL);
        reportCoalescenceResult("coalescenceResult", seqResult, f);
    }
    fprintf(f, "</coalescenceTest>\n");

    stHash_destructIterator(seqIt);
}

void compareMAFCoalescences(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash, bool onlyLeaves) {
    // Sample coalescences from the MAF (a pair of sequences from a
    // block and what genome they coalesce in)
    st_logInfo("Sampling coalescences\n");
    stSortedSet *coalescences = stSortedSet_construct3((int (*)(const void *, const void *)) coalescence_cmp, (void (*)(void *)) coalescence_destruct);
    double acceptProbability = ((double) opts->numSamples) / countPairsInMaf(opts->mafFile1, legitSequences);
    sampleCoalescences(opts->mafFile1, coalescences, acceptProbability, legitSequences, sequenceLengthHash, onlyLeaves);
    st_logInfo("Sampled %" PRIi64 " coalescences\n", stSortedSet_size(coalescences));

    st_logInfo("Finding matching coalescences\n");
    stSortedSet *matchingCoalescences = findMatchingCoalescences(opts->mafFile2, coalescences, legitSequences, onlyLeaves);
    st_logInfo("Got %" PRIi64 " matching pairs\n", stSortedSet_size(matchingCoalescences));

    st_logInfo("Accumulating results\n");
    // Overall results.
    CoalResult *aggregateResults = coalResult_init("aggregate");
    // The per-sequence result hash.
    stHash *seqResults = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void (*)(void *)) coalResult_destruct);
    buildCoalescenceResults(coalescences, matchingCoalescences, opts->speciesTree, aggregateResults, seqResults);

    FILE *outFile;
    if (opts->outFile == NULL) {
        outFile = stdout;
    } else {
        outFile = fopen(opts->outFile, "w");
    }
    reportCoalescenceResults(aggregateResults, seqResults, opts->mafFile1, opts->mafFile2, outFile);

    // Clean up.
    stSortedSet_destruct(coalescences);
    stSortedSet_destruct(matchingCoalescences);
    coalResult_destruct(aggregateResults);
    stHash_destruct(seqResults);
}
