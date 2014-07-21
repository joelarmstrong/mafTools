#ifndef __MAFPHYLOCOMPARATOR_H_
#define __MAFPHYLOCOMPARATOR_H_
#include "sonLib.h"

typedef struct {
    char *mafFile1;
    char *mafFile2;
    char *outFile;
    int64_t numSamples;
    stTree *speciesTree;
    bool onlyLeaves; // Whether the mafs have block entries for just
                     // the leaves or for ancestors as well.
} PhyloOptions;

#endif // __MAFPHYLOCOMPARATOR_H_
