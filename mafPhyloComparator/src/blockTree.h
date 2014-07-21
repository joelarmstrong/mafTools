#ifndef __BLOCKTREE_H_
#define __BLOCKTREE_H_
#include "sonLib.h"
#include "sharedMaf.h"

// For mapping between sequence headers and block tree entries.
// FIXME: don't like the name
typedef struct {
    char *species; // Sequence header, but want to use the same naming
                   // convention as the rest of mafTools.
    int64_t start; // inclusive
    int64_t end; // exclusive
    stTree *node; // node in the block tree corresponding to this row
} BlockRow;

char *parseTreeFromBlockStart(mafLine_t *line);
stTree *getMRCA(stTree *node1, stTree *node2);
stTree *getNodeFromPosition(stHash *seqToBlockRows, const char *seq, uint64_t pos);
int BlockRow_cmp(const BlockRow *row1, const BlockRow *row2);
void BlockRow_destruct(BlockRow *row);
void fillListByReversePostOrder(stTree *tree, stList *list, bool onlyLeaves);
stHash *getSeqToBlockRows(mafBlock_t *block, stTree *tree, bool onlyLeaves);
stHash *buildNameToNodeHash(stTree *tree);
bool isAncestor(char *name1, char *name2, stHash *nameToNode);

#endif
