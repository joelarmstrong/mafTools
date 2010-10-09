#include "malnJoinSets.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

/* object used to store current join block/component, which changes when new blocks are
 * merged. */
struct joinBlkComp {
    struct malnBlk *blk;
    struct malnComp *comp;
};

/* join two components updating joining object */
static void joinCompWithComp(struct joinBlkComp *joining, struct malnComp *comp2) {
    struct malnComp *joinedComp = NULL;
    struct malnBlk *joinedBlk = malnJoinBlks(joining->comp, comp2, &joinedComp);
    
    // if joining block is an a set, then just mark as done, otherwise free it.
    if (joining->blk->malnSet != NULL) {
        joining->blk->done = true;
    } else {
        malnBlk_destruct(joining->blk);
    }
    comp2->blk->done = true;
    // return resulting block for more additions
    joining->blk = joinedBlk;
    joining->comp = joinedComp;
}

/* join for specified joinable component, returning the potentially new joinedBlk.  Return
 * true if any joined. */
static bool joinCompWithSet(struct joinBlkComp *joining, struct malnSet *malnSet2) {
    bool joinedOne = false;
    stList *overComps2 = malnSet_getOverlappingPendingComps(malnSet2, joining->comp->seq, joining->comp->chromStart, joining->comp->chromEnd, mafTreeLocRoot|mafTreeLocLeaf);
    for (int i = 0; i < stList_length(overComps2); i++) {
        struct malnComp *comp2 = stList_get(overComps2, i);
        if ((!comp2->blk->done) && malnComp_canJoin(joining->comp, comp2)) {
            joinCompWithComp(joining, comp2);
            joinedOne = true;
        }
    }
    stList_destruct(overComps2);
    return joinedOne;
}

/* join blocks overlapping a block of another set. If none are joined,
 * add blk1, copy blk1 to new set */
static void joinBlkWithSet(struct malnSet *malnSetJoined, struct malnBlk *blk1, struct malnSet *malnSet2) {
    struct joinBlkComp joining = {blk1, NULL};
    // since joining creates a new block and we want to continue to search other components,
    // we start over with the first component when we create a new block.  We do this until
    // we go through all of the joinable components without actually merging a block.
    bool joinedOne = false;
    do {
        joinedOne = false;
        for (joining.comp = joining.blk->comps; (joining.comp != NULL) && (!joinedOne); joining.comp = joining.comp->next) {
            if (malnComp_joinable(joining.comp)) {
                if (joinCompWithSet(&joining, malnSet2)) {
                    joinedOne = true;
                }
            }
        }
    } while (joinedOne);
        
    if (joining.blk->malnSet == NULL) {
        // new block, something was joined
        malnSet_addBlk(malnSetJoined, joining.blk);
    } else {
        // nothing joined, copy it.
        assert(joining.blk == blk1);
        malnSet_addBlk(malnSetJoined, malnBlk_constructClone(joining.blk));
    }
}

/* add blocks that were not joined into the alignment */
static void addUndone(struct malnSet *malnSetJoined, struct malnSet *malnSet) {
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *blk;
    while ((blk = stSortedSet_getNext(iter)) != NULL) {
        if (!blk->done) {
            malnSet_addBlk(malnSetJoined, malnBlk_constructClone(blk));
            blk->done = true;
        }
    }
    stSortedSet_destructIterator(iter);
}

/* join two sets, generating a third */
struct malnSet *malnJoinSets(struct Genome *refGenome, struct malnSet *malnSet1, struct malnSet *malnSet2) {
    struct malnSet *malnSetJoined = malnSet_construct(malnSet_getGenomes(malnSet1));

    stSortedSetIterator *iter1 = malnSet_getBlocks(malnSet1);
    struct malnBlk *blk1;
    while ((blk1 = stSortedSet_getNext(iter1)) != NULL) {
        joinBlkWithSet(malnSetJoined, blk1, malnSet2);
    }
    stSortedSet_destructIterator(iter1);
    addUndone(malnSetJoined, malnSet2);
    malnSet_clearDone(malnSet1);
    malnSet_clearDone(malnSet2);
    malnSet_clearDone(malnSetJoined);
    malnSet_assert(malnSetJoined);
    return malnSetJoined;
}
