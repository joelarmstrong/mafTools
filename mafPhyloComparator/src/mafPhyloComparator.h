#ifndef __MAFPHYLOCOMPARATOR_H_
#define __MAFPHYLOCOMPARATOR_H_
typedef struct {
    char *mafFile1;
    char *mafFile2;
    int64_t numSamples;
} PhyloOptions;

// Represents an aligned pair and where their MRCA is in the gene
// tree.
typedef struct {
    char *seq1;
    int64_t pos1;
    char *seq2;
    int64_t pos2;
    char *mrca;
} Coalescence;

// FIXME: don't like the name
typedef struct {
    char *seq;
    int64_t start; // inclusive
    int64_t end; // exclusive
    stTree *node; // node in the block tree corresponding to this row
} BlockRow;
#endif // __MAFPHYLOCOMPARATOR_H_
