#include "malnJoinDups.h"
#include "malnSet.h"
#include "malnBlk.h"
#include "malnComp.h"
#include "malnDeleteBlks.h"
#include "malnJoinBlks.h"
#include "sonLibSortedSet.h"
#include "sonLibList.h"
#include "genome.h"
#include "common.h"
#include <stdbool.h>
#include <unistd.h>

// FIXME: tmp
static const bool debug = false;

/* join two blocks associated components, create and third block. Return that
 * blocks root component. Joined blocks are inserted into delete table */
static struct malnComp *joinCompWithDup(struct malnSet *malnSet, struct malnComp *comp1, struct malnComp *comp2, struct malnDeleteBlks *delBlks) {
    if (debug) {
        malnComp_dump(comp2, "joinCompWithDup comp2", stderr);
    }
    struct malnBlk *joinedBlk = malnJoinBlks(comp1, comp2, NULL);
    malnDeleteBlks_flag(delBlks, comp1->blk);
    malnDeleteBlks_flag(delBlks, comp2->blk);
    return malnBlk_getRootComp(joinedBlk);
}

/* attempt to join a block with other blocks using the specified component.
 * Return updated block when one join is achieved, or NULL if no join was done. */
static struct malnBlk *joinCompWithDups(struct malnSet *malnSet, struct malnBlk *joinBlk, struct malnComp *joinComp, struct malnDeleteBlks *delBlks) {
    if (debug) {
        malnComp_dump(joinComp, "joinCompWithDups", stderr);
    }
    joinComp->blk->done = true;
    stList *overComps = malnSet_getOverlappingPendingComps(malnSet, joinComp->seq, joinComp->chromStart, joinComp->chromEnd, mafTreeLocAll);
    for (int i = 0; i < stList_length(overComps); i++) {
        struct malnComp *dupComp = stList_get(overComps, i);
        if (!malnDeleteBlks_contains(delBlks, dupComp->blk)) {
            joinComp = joinCompWithDup(malnSet, joinComp, dupComp, delBlks);
            joinComp->blk->done = true;
        }
    }
    stList_destruct(overComps);
    return (joinComp->blk != joinBlk) ? joinComp->blk : NULL;
}

/* Join one block with any duplications of that block in the same
 * set. Duplicates are added to a table and skipped so that iterators are not
 * invalidated. */
static void joinBlkWithDups(struct malnSet *malnSet, struct malnBlk *joinBlk, struct malnDeleteBlks *delBlks) {
    // iterate over each component in joinBlk, looking for overlapping components
    // in other blocks.  A join creates a new block, so we start over until no
    // block is joined with this block.

    if (debug) { // FIXME: tmp
        malnBlk_dump(joinBlk, "joinBlkWithDups", stderr);
    }
    bool joinedSome = FALSE;
    struct malnBlk *newJoinBlk;
    do {
        newJoinBlk = NULL;
        for (struct malnComp *joinComp = joinBlk->comps; (joinComp != NULL) && (newJoinBlk == NULL); joinComp = joinComp->next) {
            newJoinBlk = joinCompWithDups(malnSet, joinBlk, joinComp, delBlks);
        }
        if (newJoinBlk != NULL) {
            joinBlk = newJoinBlk;
            joinedSome = true;
        }
    } while (newJoinBlk != NULL);

    if (joinedSome) {
        malnSet_addBlk(malnSet, joinBlk);
    }
}

/* Join duplication blocks in a set, which evolver outputs as separate
 * blocks. */
void malnJoin_joinSetDups(struct malnSet *malnSet) {
    struct malnDeleteBlks *delBlks = malnDeleteBlks_construct();
    stSortedSetIterator *iter = malnSet_getBlocks(malnSet);
    struct malnBlk *joinBlk;
    while ((joinBlk = stSortedSet_getNext(iter)) != NULL) {
        if (!malnDeleteBlks_contains(delBlks, joinBlk)) {
            joinBlkWithDups(malnSet, joinBlk, delBlks);
        }
    }
    stSortedSet_destructIterator(iter);
    malnDeleteBlks_destruct(delBlks);
    malnSet_assert(malnSet);
    malnSet_clearDone(malnSet);
}

