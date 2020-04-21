//
// Created by alberto on 19/02/2020.
//

#ifndef TBKP_TBKP_UB_H
#define TBKP_TBKP_UB_H

#include "tbkp_instance.h"
#include "tbkp_solution.h"
#include "tbkp_bb_fixed_status.h"
#include "tbkp_bb_stats.h"
#include "tbkp_bb_defs.h"
#include <stddef.h>

/** Solves the Time-Bomb Knapsack Problem by branch and bound.
 *
 * @param instance  The TBKP instance we are solving.
 * @param stats     Object containing parameters, and where to save solution statistics.
 * @return          A pointer to the solution of the problem.
 */
TBKPSolution* tbkp_branch_and_bound(const TBKPInstance* instance, TBKPBBStats* stats);

#endif //TBKP_TBKP_UB_H
