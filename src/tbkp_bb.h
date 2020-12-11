//
// Created by alberto on 19/02/2020.
//

#ifndef TBKP_TBKP_BB_H
#define TBKP_TBKP_BB_H

#include "tbkp_instance.h"
#include "tbkp_bb_solution.h"
#include "tbkp_bb_fixed_status.h"
#include "tbkp_bb_stats.h"
#include "tbkp_bb_defs.h"
#include <stddef.h>
#include <gurobi_c.h>

typedef struct {
    const TBKPInstance* instance;
    const TBKPParams* params;
    TBKPBBSolution* solution;
    TBKPBBFixedStatus* x;
    TBKPBBStats* stats;
    size_t* n_nodes;
    GRBmodel* boole_grb_model;
} TBKPBBAlgStatus;

typedef struct {
    double prod_probabilities;
    uint_fast32_t sum_profits;
    uint_fast32_t res_capacity;
} TBKPBBResidualInstance;

/** Solves the Time-Bomb Knapsack Problem by branch and bound.
 *
 * @param instance  The TBKP instance we are solving.
 * @param stats     Out-param where to store solution statistics.
 * @param params    Struct containing user-set parameters.
 * @return          A pointer to the solution of the problem.
 */
TBKPBBSolution* tbkp_branch_and_bound(const TBKPInstance* instance, TBKPBBStats* stats, TBKPParams* params);

double solve_cont(const TBKPInstance *const instance, TBKPBBFixedStatus* xbra, double *xc);


#endif //TBKP_TBKP_BB_H
