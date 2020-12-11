//
// Created by alberto on 23/02/2020.
//

#ifndef TBKP_TBKP_BOOLE_SOL_H
#define TBKP_TBKP_BOOLE_SOL_H

#include "tbkp_instance.h"
#include "tbkp_params.h"
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
     *  any solution can give a LB.
     *
     *  The objective function of the Boole BQ problem is a lower bound. However,
     *  we can obtain a tighter LB by computing the original TBKP objective
     *  function using the items packed in the Boole BQ problem, which is what
     *  we store in this variable.
     */
    double lb;

    /** Part of the LB coming from the sum of the profits. */
    uint_fast32_t lb_sum_profits;

    /** Part of the LB coming from the product of the probabilities. */
    double lb_product_probabilities;

    /** Number of items packed in the solution of the Boole BQ problem. */
    size_t n_items;

    /** Items packed in the solution of the Boole BQ problem. */
    size_t* items;

    /** Time in seconds to compute the solution. */
    float time_to_compute;
} TBKPBooleSol;

/** Generate a base model for solving the Boole BQ problem as a quadratic problem.
 * @param instance Const pointer to a Time-bomb Knapsack problem.
 * @return         A Gurobi model to be reused for the solution of Boole BQ problems.
 */
GRBmodel* tbkp_boolesol_quad_base_model(const TBKPInstance* instance);

/** Gets the lower bound from solving the Boole BQ problem as-is, i.e., as
 *  a quadratic problem using the Gurobi optimiser.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack problem.
 * @param params    Solver parameters.
 * @param grb_model Base Gurobi model.
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
TBKPBooleSol tbkp_boolesol_quad_gurobi_get(
        const TBKPInstance* instance,
        const TBKPParams* params,
        GRBmodel* grb_model,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity);

/** Generate a base model for solving the Boole BQ problem using a linearisation.
 * @param instance Const pointer to a Time-bomb Knapsack problem.
 * @return         A Gurobi model to be reused for the solution of Boole BQ problems.
 */
GRBmodel* tbkp_boolesol_lin_base_model(const TBKPInstance* instance);

/** Gets the lower bound from solving the Boole BQ problem using a linearisation
 *  and getting its solution via the Gurobi solver.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack problem.
 * @param params    Solver parameters.
 * @param grb_model Base Gurobi model.
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
TBKPBooleSol tbkp_boolesol_lin_gurobi_get(
        const TBKPInstance* instance,
        const TBKPParams* params,
        GRBmodel* grb_model,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity);

/** Computes the original TBKP objective value of a solution obtained
 *  using the Boole BQ problem.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack problem.
 * @param sol       Boole BQ solution.
 */
void tbkp_boolesol_compute_exact_obj(
        const TBKPInstance* instance,
        TBKPBooleSol* sol);

void tbkp_boolesol_print(const TBKPBooleSol* sol);
void tbkp_boolesol_free_inside(TBKPBooleSol* boolesol_ptr);

#endif //TBKP_TBKP_BOOLE_SOL_H
