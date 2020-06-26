#ifndef TBKP_TBKP_CR_SOL_H
#define TBKP_TBKP_CR_SOL_H

#include <stddef.h>
#include "tbkp_instance.h"
#include "tbkp_bb_fixed_status.h"

/**
 * Structure modelling the result of solving the continuous relaxation of the TB-01-KP.
 */
typedef struct {
    /** Upper bound, given by the optimal solution of the CR. */
    float ub;

    /** If the solution of the continuous relaxation is integer, the UB is also a LB.
     */
    _Bool is_lb;

    /** If the solution is integer, part of the LB coming from summing the profits. */
    uint_fast32_t lb_sum_profits;

    /** If the solution is integer, part of the LB coming from multiplying the probabilities. */
    float lb_product_probabilities;

    /** If the solution is integer, this is the number of items packed. */
    size_t n_items;

    /** If the solution is integer, these are the items packed.
     */
    size_t* items;

    /** Time in seconds to compute the solution. */
    float time_to_compute;
} TBKPContinuousRelaxationSol;

/** Gets the continuous relaxation bound, solving a convex optimisation problem.
 *
 * @param instance  Const pointer to a Time-bomb Knapsack instance.
 * @param xbra      Fixed status of the variables (fixed in the knapsack; excluded; unfixed).
 * @return          The bound coming from the continuous relaxation.
 */
TBKPContinuousRelaxationSol tbkp_crsol_get(
        const TBKPInstance* instance,
        const TBKPBBFixedStatus* xbra);

void tbkp_crsol_free_inside(TBKPContinuousRelaxationSol* desol_ptr);
void tbkp_crsol_print(const TBKPContinuousRelaxationSol* sol);

#endif //TBKP_TBKP_DE_SOL_H
