#ifndef __COALESCENCES_H_
#define __COALESCENCES_H_
#include "sonLib.h"
#include "sharedMaf.h"
#include "mafPhyloComparator.h"

// Represents an aligned pair and where their MRCA is in the gene
// tree.
typedef struct {
    char *seq1;
    int64_t pos1;
    char *seq2;
    int64_t pos2;
    char *mrca;
} Coalescence;

// result structure for coalescence result.
typedef struct {
    char *seq; // == "aggregate" for overall results
    int64_t sampledPairs;
    uint64_t matchingCoalescences; // # of coalescences whose (seq,
                                   // pos) pairs match
    uint64_t identicalCoalescences; // # of matching coalescences that
                                    // # coalesce at the same place
    uint64_t earlyCoalescences; // Coalesced early in B relative to A.
    uint64_t lateCoalescences; // Coalesced late in B relative to A.
} CoalResult;

int coalescence_cmp(const Coalescence *coal1, const Coalescence *coal2);
void coalescence_destruct(Coalescence *coal);
void coalescencesFromPairs(stSortedSet *pairs, stHash *seqToBlockRows, stSortedSet *coalescences);
void walkBlockSamplingCoalescences(char *mafFileName, mafBlock_t *block, stSortedSet *coalescences, double acceptProbability, stSet *legitSequences, stHash *sequenceLengthHash, uint64_t *chooseTwoArray, bool onlyLeaves);

// Sample, compare and report coalescences from two MAFs.
void compareMAFCoalescences(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash, bool onlyLeaves);
#endif
