//
// Created by alberto on 19/02/2020.
//

#ifndef TBKP_TBKP_UB_H
#define TBKP_TBKP_UB_H

#include "tbkp_instance.h"
#include <stddef.h>
#include <stdbool.h>
#include <time.h>

#define NO_TIMEBOMB_ITEM_TO_BRANCH -1

/** Represents a solution to the Time-Bomb Knapsack Problem. */
typedef struct {
    /** Vector of bools. x[i] is true iff item i is packed. */
    _Bool* x;

    /** Product-of-probabilities component of the objective function. */
    float prod_probabilities;

    /** Sum-of-profits component of the objective function. */
    uint_fast32_t sum_profits;

    /** Objective value. */
    float value;
} TBKPSolution;

TBKPSolution* tbkp_sol_init(const TBKPInstance* instance);
void tbkp_sol_print(const TBKPSolution *const solution, const TBKPInstance* instance);
void tbkp_sol_free(TBKPSolution** solution_ptr);

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

/** Collects statistics on the B&B solution process. */
typedef struct {
    clock_t start_time;
    clock_t end_time;
    float timeout;
    float elapsed_time;
    float ub; /* Starts at value -1.0f */
    float lb;
    float gap;
    size_t n_nodes;
} TBKPStats;

/** Initialises an empty statistics object. */
TBKPStats tbkp_stats_init(float timeout);

/** Prints a summary of the statistics/ */
void tbkp_stats_print(const TBKPStats *const stats);

/** Solves the Time-Bomb Knapsack Problem by branch and bound.
 *
 * @param instance  The TBKP instance we are solving.
 * @param stats     Object where to save solution statistics.
 * @return          A pointer to the solution of the problem.
 */
TBKPSolution* tbkp_branch_and_bound(const TBKPInstance* instance, TBKPStats* stats);

/** Finds the next time-bomb item to branch on. If no such item exists (because we either
 *  ran out of TB items, or there is no item which fits in the residual capacity), returns
 *  NO_TIMEBOMB_ITEM_TO_BRANCH.
 *
 * @param instance  The TBKP instance we are solving.
 * @param x         An array keeping track of the branching decisions:
 *                  x[i] == FIXED_DONT_PACK if item i has been fixed to NOT being packed
 *                  x[i] == FIXED_PACK if item i has been fixed to being packed
 *                  x[i] == UNFIXED if item i is unfixed.
 * @param capacity  Residual capacity of the knapsack.
 * @return          The index of the item to branch on, or NO_TIMEBOMB_ITEM_TO_BRANCH.
 */
int tbkp_bb_branch_item(const TBKPInstance* instance, const TBKPBBFixedStatus* x, uint_fast32_t capacity, uint_fast32_t sum_profits);

/** Solves a Branch-and-Bound node and possible updates the current solution.
 *
 * @param instance              TBKP instance we are solving.
 * @param nnodes                Number of nodes in the BB tree.
 * @param x                     Vector keeping track of fixed/unfixed tb items.
 * @param prod_probabilities    Product-of-probabilities part of the objective function.
 * @param sum_profits           Sum-of-profits part of the objective function.
 * @param res_capacity          Residual capacity of the knapsack at this node.
 * @param solution              Current solution.
 */
void tbkp_bb_solve_node(
        const TBKPInstance* instance,
        size_t* nnodes,
        TBKPBBFixedStatus* x,
        float parent_ub,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity,
        TBKPSolution *solution,
        TBKPStats* stats);

#endif //TBKP_TBKP_UB_H
