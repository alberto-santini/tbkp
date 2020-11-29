//
// Created by alberto on 21/04/2020.
//

#ifndef TBKP_TBKP_BB_SOLUTION_H
#define TBKP_TBKP_BB_SOLUTION_H

#include "tbkp_instance.h"
#include "tbkp_bb_stats.h"
#include <stdint.h>

/** Represents a solution to the Time-Bomb Knapsack Problem. */
typedef struct {
    /** Vector of bools. x[i] is true iff item i is packed. */
    _Bool* x;

    /** Product-of-probabilities component of the objective function. */
    double prod_probabilities;

    /** Sum-of-profits component of the objective function. */
    uint_fast32_t sum_profits;

    /** Objective value. */
    double value;
} TBKPBBSolution;

TBKPBBSolution* tbkp_sol_init(const TBKPInstance* instance);
void tbkp_sol_print(const TBKPBBSolution* solution, const TBKPInstance* instance);
void tbkp_sol_free(TBKPBBSolution** solution_ptr);

#endif //TBKP_TBKP_BB_SOLUTION_H
