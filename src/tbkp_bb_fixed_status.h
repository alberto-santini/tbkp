//
// Created by alberto on 21/04/2020.
//

#ifndef TBKP_TBKP_BB_FIXED_STATUS_H
#define TBKP_TBKP_BB_FIXED_STATUS_H

/** Represents the current status of an item in the branch-and-bound tree.
 * If UNFIXED, then the item is free (can be packed or not).
 * If FIXED_DONT_PACK, then the item is excluded from the knapsack.
 * If FIXED_PACK, then the item is forced into the knapsack.
 * Only applies to time-bomb items.
 */
typedef enum {
    UNFIXED = -1,
    FIXED_DONT_PACK = 0,
    FIXED_PACK = 1
} TBKPBBFixedStatus;

#endif //TBKP_TBKP_BB_FIXED_STATUS_H
