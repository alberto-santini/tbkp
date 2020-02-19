//
// Created by alberto on 19/02/2020.
//

#ifndef TBKP_TBKP_UB_H
#define TBKP_TBKP_UB_H

#include <stddef.h>
#include "tbkp_instance.h"

/**
 * Structure modelling the result of solving the "deterministic equivalent" 0-1 Knapsack.
 * This is the problem where the profit of each object is replaced by profit*probability;
 * then the problem is solved like a deterministic 0-1 Knapsack problem.
 */
typedef struct {
    /** Upper bound, given by the optimal solution of the deterministic equivalent 0-1KP. */
    float ub;

    /** Number of items packed in the solution of the deterministic equivalent 0-1 KP. */
    size_t n_items;

    /** Items in the optimal solution of the deterministic equivalent 0-1KP.
     *  They can also be used to get a feasible solution to the TBKP.
     */
    size_t* items;
} TBKPDeterministicEqUB;

/** Gets the "deterministic equivalent" upper bound, solving a 0-1 Knapsack with COMBO.
 *
 * @param instance          Const pointer to a Time-bomb Knapsack instance.
 * @param n_items           Number of items of the instance that we consider. This allows to get an upper bound only
 *                          on a "sub-instance" in which some items have been discarded.
 * @param items             Vector of items that we want to consider.
 * @param capacity          Capacity of the deterministic knapsack. This parameter, together with n_items and items
 *                          allows to solve subproblems in which some original item is "fixed" in the knapsack.
 * @return                  The value of the UB and the objects packed in the deterministic 0-1 Knapsack used.
 */
TBKPDeterministicEqUB tbkp_deub_get(
        const TBKPInstance* instance,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity);
void tbkp_deub_free_inside(TBKPDeterministicEqUB* deub_ptr);
void tbkp_deub_print(const TBKPDeterministicEqUB* ub);

#endif //TBKP_TBKP_UB_H
