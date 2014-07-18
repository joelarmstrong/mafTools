#ifndef __MAFPHYLOCOMPARATOR_H_
#define __MAFPHYLOCOMPARATOR_H_
typedef struct {
    char *mafFile1;
    char *mafFile2;
    char *outFile;
    int64_t numSamples;
    stTree *speciesTree;
} PhyloOptions;

#endif // __MAFPHYLOCOMPARATOR_H_
