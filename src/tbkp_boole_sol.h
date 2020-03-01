//
// Created by alberto on 23/02/2020.
//

#ifndef TBKP_TBKP_BOOLE_SOL_H
#define TBKP_TBKP_BOOLE_SOL_H

#include "tbkp_instance.h"
#include <stddef.h>
#include <gurobi_c.h>

extern GRBenv* grb_env;

/**
 * Structure modelling the solution to the Binary Quadratic problem obtained
 * using Boole's law on the original 0-1 Time-bomb Knapsack objective function.
 */
typedef struct {
    /** Lower bound given by any solution to the Boole BQ problem. The optimal
     *  solution gives the tightest bound, but it's not necessary to have it:
     *  any solution gives a LB.
     *
     *  The objective function of the Boole BQ problem is a lower bound. Additionally,
     *  one can compute the 01-TB-KP objective value on the items packed by
     *  the Boole BQ solution. At the end, one can take the highest of the two.
     */
    float lb;

    /** Number of items packed in the solution of the Boole BQ problem. */
    size_t n_items;

    /** Items packed in the solution of the Boole BQ problem. */
    size_t* items;
} TBKPBoolSol;

/** Gets the lower bound from solving the Boole BQ problem as-is, i.e., as
 *  a quadratic problem using the Gurobi optimiser.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack problem.
 * @param n_items   Number of items to consider when solving the Boole BQ
 *                  problem. This can be smaller than the number of items in
 *                  the TBKP instance, if some items are excluded or fixed,
 *                  e.g., because we are exploring a B&B tree.
 * @param items     List of items from the original TBKP instance, which we are
 *                  considering when solving the Boole BQ problem.
 * @param capacity  Capacity of the Boole BQ knapsack (which can be smaller
 *                  than the original TBKP instance's capacity).
 * @return          The value of the LB and the objects selected by the Boole
 *                  QB problem.
 */
TBKPBoolSol tbkp_boolsol_quad_gurobi_get(
        const TBKPInstance* instance,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity);

/** Gets the lower bound from solving the Bool BQ problem using a linearisation
 *  and getting its solution via the Gurobi solver.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack problem.
 * @param n_items   Number of items to consider when solving the Boole BQ
 *                  problem. This can be smaller than the number of items in
 *                  the TBKP instance, if some items are excluded or fixed,
 *                  e.g., because we are exploring a B&B tree.
 * @param items     List of items from the original TBKP instance, which we are
 *                  considering when solving the Boole BQ problem.
 * @param capacity  Capacity of the Boole BQ knapsack (which can be smaller
 *                  than the original TBKP instance's capacity).
 * @return          The value of the LB and the objects selected by the Boole
 *                  QB problem.
 */
TBKPBoolSol tbkp_boolsol_lin_gurobi_get(
        const TBKPInstance* instance,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity);

/** Computes the original TBKP objective value of a solution obtained
 *  using the Boole BQ problem.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack problem.
 * @param sol       Boole BQ solution.
 */
void tbkp_boolsol_compute_exact_obj(
        const TBKPInstance* instance,
        TBKPBoolSol* sol);
void tbkp_boolsol_print(const TBKPBoolSol* sol);
void tbkp_boolsol_free_inside(TBKPBoolSol* boolsol_ptr);

#endif //TBKP_TBKP_BOOLE_SOL_H
