#ifndef __COALESCENCES_H_
#define __COALESCENCES_H_
#include "sonLib.h"
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

// Sample, compare and report coalescences from two MAFs.
void compareMAFCoalescences(PhyloOptions *opts, stSet *legitSequences, stHash *sequenceLengthHash);
#endif
