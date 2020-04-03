//
// Created by alberto on 19/02/2020.
//

#ifndef TBKP_TBKP_DE_SOL_H
#define TBKP_TBKP_DE_SOL_H

#include <stddef.h>
#include "tbkp_instance.h"

/**
 * Structure modelling the result of solving the "deterministic equivalent" 0-1 Knapsack.
 * This is the problem where the profit of each object is replaced by profit*probability;
 * then the problem is solved like a deterministic 0-1 Knapsack problem.
 */
typedef struct {
    /** Upper bound, given by the optimal solution of the deterministic equivalent 0-1 KP. */
    float ub;

    /** Lower bound, obtained calculating the TBKP objective function on the solution
     *  returned by the deterministic equivalent 0-1 KP.
     */
    float lb;

    /** Part of the LB coming from summing the profits. */
    uint_fast32_t lb_sum_profits;

    /** Part of the LB coming from multiplying the probabilities. */
    float lb_product_probabilities;

    /** Number of items packed in the solution of the deterministic equivalent 0-1 KP. */
    size_t n_items;

    /** Items in the optimal solution of the deterministic equivalent 0-1KP.
     *  They can also be used to get a feasible solution to the TBKP.
     */
    size_t* items;
} TBKPDeterministicEqSol;

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
TBKPDeterministicEqSol tbkp_desol_get(
        const TBKPInstance* instance,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity);

void tbkp_desol_free_inside(TBKPDeterministicEqSol* desol_ptr);
void tbkp_desol_print(const TBKPDeterministicEqSol* sol);

#endif //TBKP_TBKP_DE_SOL_H
