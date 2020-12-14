//
// last update: Michele 03/03/2020.
//

#include "tbkp_bb.h"
#include "tbkp_de_sol.h"
#include "tbkp_cr_sol.h"
#include "tbkp_boole_sol.h"
#include "utils/pdqsort_c.h"
#include <stdlib.h>
#include <stdio.h>
#include <combo.h>
#include <assert.h>
#include <float.h>
#include <time.h>

#ifndef BB_VERBOSITY_CURRENT
#define BB_VERBOSITY_CURRENT 1001
#endif

#define BB_VERBOSITY_INFO 1000
#define EPS 1e-6

/************************************************************
 * LOCAL HELPER FUNCTIONS                                   *
 ************************************************************/

static void tbkp_bb_solve_node(
        TBKPBBAlgStatus* status, double parent_ub, TBKPBBResidualInstance residual, _Bool early_combo);

/** Finds the next time-bomb item to branch on. If no such item exists (because we either
 *  ran out of TB items, or there is no item which fits in the residual capacity), returns
 *  NO_TIMEBOMB_ITEM_TO_BRANCH.
 *
 * @param status      Current state of the algorithm.
 * @param residual    Current residual instance.
 * @return            The index of the item to branch on, or NO_TIMEBOMB_ITEM_TO_BRANCH.
 */
static int tbkp_bb_branch_item(const TBKPBBAlgStatus *const status, const TBKPBBResidualInstance *const residual) {
    // Find the first uncertain item that is not fixed and fits in the residual capacity
    for(size_t i = 0u; i < status->instance->n_items; ++i) {
        if(status->instance->probabilities[i] >= 1.0 - EPS) continue;
        if(status->x[i] != UNFIXED) continue;
        if(status->instance->weights[i] > residual->res_capacity) continue;

        if(status->params->use_early_pruning) {
            // Pruning: skip branch in case the solution value cannot increase
            double scorej = ((double)status->instance->profits[i] * status->instance->probabilities[i]) /
                    (1.0 - status->instance->probabilities[i]);

            if(scorej < residual->sum_profits) {
                continue;
            }
        }

        return (int)i;
    }

    return NO_TIMEBOMB_ITEM_TO_BRANCH;
}

/**
 * Checks if timeout occurred while exploring the B&B tree.
 * It updates the TBKPBBStats object with
 * @param status        Current state of the algorithm.
 * @param parent_ub     Upper bound at the parent node.
 * @param force_timeout Forces assuming a timeout even if not true. In this way, we can
 *                      recycle this function for the case in which we reached the max
 *                      number of nodes to explore, as well as for real timeouts.
 * @return              Returns true iff a timeout occurred.
 */
static _Bool timeout(TBKPBBAlgStatus* status, double parent_ub, _Bool force_timeout) {
    clock_t current_time = clock();
    float el_time = (float)(current_time - status->stats->start_time) / CLOCKS_PER_SEC;

    if(force_timeout || el_time > status->params->timeout) {
        // Upon timeout, the UB is the worst UB of the open nodes.
        if( parent_ub != INITIAL_UB_PLACEHOLDER &&
            (status->stats->ub == INITIAL_UB_PLACEHOLDER || parent_ub > status->stats->ub + EPS))
        {
            status->stats->ub = parent_ub;
        }

        return true;
    }

    return false;
}

/**
 * Returns the number of unfixed TB items at current node.
 * @param instance  Current state of the algorithm.
 * @return          The number of items with UNFIXED status.
 */
static size_t num_unfixed_items(const TBKPBBAlgStatus *const status) {
    size_t n_unfixed_items = 0u;
    for(size_t j = 0u; j < status->instance->n_items; ++j) {
        if (status->x[j] == UNFIXED) {
            n_unfixed_items++;
        }
    }
    return n_unfixed_items;
}

/**
 * Gives the residual instance at the current node, i.e., the instance without
 * all fixed TB items.
 *
 * @param status            Current state of the algorithm.
 * @param n_unfixed_items   Number of unfixed items in list x.
 * @return                  An array with the indices of items in the residual instance.
 */
static size_t* residual_instance(const TBKPBBAlgStatus *const status, size_t n_unfixed_items) {
    size_t* items = malloc(n_unfixed_items * sizeof(*items));
    if(!items) {
        printf("Cannot allocate memory for the residual instance's items\n");
        exit(EXIT_FAILURE);
    }

    size_t cnt = 0;
    for(size_t j = 0; j < status->instance->n_items; j++) {
        if(status->x[j] == UNFIXED) {
            items[cnt++] = j;
        }
    }

    return items;
}

/**
 * Structure to contain info from the Deterministic Equivalent bounds,
 * which are useful in the B&B nodes. It contains the UB and LB at the
 * local (current) node, plus a flag indicating if the current node should
 * be prunned - i.e., because the local UB is worse than the best-known
 * LB.
 */
typedef struct {
    double local_ub;
    double local_lb;
    _Bool should_prune;
} DEBounds;

/**
 * Structure to contain infor from the Continuous Relaxation bound,
 * which are useful in the B&B nodes. It contains the UB at the local
 * node, a flag telling whether the UB is also a LB (i.e., the solution
 * was integer), and a flag indicating if the current node should be
 * pruned - i.e., because the local UB is worse than the best-known LB.
 */
typedef struct {
    double local_ub;
    _Bool is_lb;
    _Bool should_prune;
} CRBound;

/**
 * Updates the current best solution from the LB-solution gotten from the
 * Deterministic Equivalent solution + the fixed items.
 *
 * @param status                Current state of the algorithm.
 * @param residual              Current residual instance.
 * @param desol                 Deterministic Equivalent solution.
 * @param local_lb              Local LB, i.e., obj value of the DE solution.
 */
static void update_best_solution_from_de(
        TBKPBBAlgStatus* status,
        const TBKPBBResidualInstance *const residual,
        const TBKPDeterministicEqSol *const desol,
        double local_lb
) {
    status->solution->value = local_lb;
    status->solution->prod_probabilities = desol->lb_product_probabilities * residual->prod_probabilities;
    status->solution->sum_profits = desol->lb_sum_profits + residual->sum_profits;

    for(size_t i = 0u; i < status->instance->n_items; ++i) {
        if(status->x[i] != UNFIXED) {
            status->solution->x[i] = status->x[i];
        }
    }

    for(size_t i = 0u; i < desol->n_items; ++i) {
        status->solution->x[desol->items[i]] = true;
    }

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\t\tSolution update from DE: new value %f\n", status->solution->value);
    }
}

/**
 * Updates the current best solution from the LB-solution gotten from the
 * Deterministic Equivalent solution only, disregarding the fixed items.
 *
 * @param status                Current state of the algorithm.
 * @param desol                 Deterministic Equivalent solution.
 */
static void update_best_solution_from_de_alone(
        TBKPBBAlgStatus* status,
        const TBKPDeterministicEqSol *const desol
) {
    status->solution->value = desol->lb;
    status->solution->prod_probabilities = desol->lb_product_probabilities;
    status->solution->sum_profits = desol->lb_sum_profits;

    for(size_t i = 0u; i < status->instance->n_items; ++i) {
        status->solution->x[i] = false;
    }

    for(size_t i = 0u; i < desol->n_items; ++i) {
        status->solution->x[desol->items[i]] = true;
    }

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\t\tSolution update from DE alone: new value %f\n", status->solution->value);
    }
}

/**
 * Updates the current best solution from the Continuous Relaxation solution.
 *
 * @param status                Current state of the algorithm.
 * @param crsol                 Continuous Relaxation solution.
 * @param local_lb              Local LB, i.e., obj value of the DE solution.
 */
static void update_best_solution_from_cr(
        TBKPBBAlgStatus* status,
        const TBKPContinuousRelaxationSol *const crsol,
        double local_lb
) {
    status->solution->value = local_lb;
    status->solution->prod_probabilities = crsol->lb_product_probabilities;
    status->solution->sum_profits = crsol->lb_sum_profits;

    for(size_t i = 0u; i < status->instance->n_items; ++i) {
        if(status->x[i] != UNFIXED) {
            status->solution->x[i] = status->x[i];
        }
    }

    for(size_t i = 0u; i < crsol->n_items; ++i) {
        status->solution->x[crsol->items[i]] = true;
    }

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\t\tSolution update from CR: new value %f\n", status->solution->value);
    }
}

/**
 * Updates the current best solution from the LB-solution gotten from the
 * Boole bound solution.
 *
 * @param status                Current state of the algorithm.
 * @param residual              Current residual instance.
 * @param boolesol              Boole solution.
 * @param local_lb              Local LB, i.e., obj value of the Boole solution.
 */
static void update_best_solution_from_boole(
        TBKPBBAlgStatus* status,
        const TBKPBBResidualInstance *const residual,
        const TBKPBooleSol *const boolesol,
        double local_lb
) {
    status->solution->value = local_lb;
    status->solution->prod_probabilities = boolesol->lb_product_probabilities * residual->prod_probabilities;
    status->solution->sum_profits = boolesol->lb_sum_profits + residual->sum_profits;

    for(size_t i = 0u; i < status->instance->n_items; ++i) {
        if(status->x[i] != UNFIXED) {
            status->solution->x[i] = status->x[i];
        }
    }

    for(size_t i = 0u; i < boolesol->n_items; ++i) {
        status->solution->x[boolesol->items[i]] = true;
    }

    if(BB_VERBOSITY_CURRENT > BB_VERBOSITY_INFO) {
        printf("\t\tSolution update from Boole: new value %f\n", status->solution->value);
    }
}

/** Gets the Continuous Relaxation UB at the current node.
 *  If the UB is worse than the current best LB, it signals to prune the current node.
 *  If the relaxation solution is integer, the UB also works as a LB.
 *  In this case, if the LB is better than the current best LB, it updates it.
 *
 *  @param status               Current state of the algorithm.
 *  @return                     The CRBounds object containind the CR bounds and the pruning flag.
 */
static CRBound get_cr_bound(
    TBKPBBAlgStatus *const status,
    size_t current_node
) {
    // Compute the continuous relaxation bound
    TBKPContinuousRelaxationSol crsol = tbkp_crsol_get(status->instance, status->x);

    ++(status->stats->n_cr_called);
    status->stats->tot_time_cr += crsol.time_to_compute;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\t\tUB (continuous relaxation): %f (vs %f)\n", crsol.ub, status->solution->value);
    }

    if(crsol.ub <= status->solution->value) {
        return (CRBound){.local_ub = crsol.ub, .is_lb = crsol.is_lb, .should_prune = true};
    }

    if(current_node == 1u) {
        status->stats->cr_ub_at_root = crsol.ub;
        status->stats->time_to_compute_cr_at_root = crsol.time_to_compute;
    }

    // Check if the solution of the continuous relaxation is integer and improves
    if(crsol.is_lb) {
        const double local_lb = crsol.ub;

        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("\t\tLB (continuous relaxation): %f (vs %f)\n", crsol.ub, status->solution->value);
        }

        // Possibly update the best integer solution
        if(local_lb > status->stats->lb) {
            assert(local_lb > status->solution->value);

            status->stats->lb = local_lb;
            update_best_solution_from_cr(status, &crsol, local_lb);
        }
    }

    tbkp_crsol_free_inside(&crsol);

    // If the cont relax is integer (crsol.is_lb) we can prune the node.
    return (CRBound){.local_ub = crsol.ub, .is_lb = crsol.is_lb, .should_prune = crsol.is_lb};
}

/**
 * Gets the Deterministic Equivalent bounds (UB and LB) at the current node.
 * If the UB is worse than the current best LB, it signals to prune the current node.
 * It the LB is better than the current best LB, it updates it.
 *
 * @param status                Current state of the algorithm.
 * @param residual              Residual instance at the current node.
 * @param items                 List with the indices of items in the current residual instance.
 * @param n_unfixed_items       Number of unfixed items in the current residual instance.
 * @return                      The DEBounds object containing the DE bounds and the pruning flag.
 */
static DEBounds get_de_bounds(
        TBKPBBAlgStatus *const status,
        const TBKPBBResidualInstance *const residual,
        size_t current_node,
        size_t* items,
        size_t n_unfixed_items
) {
    // Solve the deterministic relaxation
    TBKPDeterministicEqSol desol = tbkp_desol_get(status->instance, n_unfixed_items, items, residual->res_capacity);

    ++(status->stats->n_de_called);
    status->stats->tot_time_de += desol.time_to_compute;

    // Compute the local upper bound and possibly kill the node
    double local_ub = ((double)residual->sum_profits + desol.ub) * residual->prod_probabilities;
    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\t\tUB (deterministic relaxation): %f\n", local_ub);
    }

    if(local_ub <= status->solution->value) {
        tbkp_desol_free_inside(&desol);
        free(items); items = NULL;
        return (DEBounds){.local_ub = local_ub, .local_lb = 0.0f, .should_prune = true};
    }
    
    double local_lb = (double) (desol.lb_sum_profits + residual->sum_profits) *
                     (desol.lb_product_probabilities * residual->prod_probabilities);

    // If the heuristic solution obtained with the deterministic bound is better
    // by itself than it is when combining it with the residual instance, we should
    // instead use this one as the best integer solution produced by this bound.
    _Bool better_de_alone = false;
    if(desol.lb > local_lb) {
        better_de_alone = true;
        local_lb = desol.lb;
    }

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\t\tLB (deterministic relaxation): %f (vs %f)\n", local_lb, status->solution->value);
    }

    if(current_node == 1u) {
        status->stats->de_lb_at_root = local_lb;
        status->stats->de_ub_at_root = local_ub;
        status->stats->time_to_compute_de_at_root = desol.time_to_compute;
    }

    // Possibly update the best integer solution
    if(local_lb > status->stats->lb) {
        assert(local_lb > status->solution->value);

        status->stats->lb = local_lb;

        if(better_de_alone) {
            update_best_solution_from_de_alone(status, &desol);
        } else {
            update_best_solution_from_de(status, residual, &desol, local_lb);
        }
    }

    tbkp_desol_free_inside(&desol);

    return (DEBounds){.local_ub = local_ub, .local_lb = local_lb, .should_prune = false};
}

/**
 * Gets the Boole bound (LB) at the current node.
 * It the LB is better than the current best LB, it updates it.
 *
 * @param status                Current state of the algorithm.
 * @param residual              Residual instance at the current node.
 * @param items                 List with the indices of items in the current residual instance.
 * @param n_unfixed_items       Number of unfixed items in the current residual instance.
 * @return                      The lower bound.
 */
double get_boole_bound(
        TBKPBBAlgStatus *const status,
        const TBKPBBResidualInstance *const residual,
        size_t current_node
) {
    TBKPBooleSol boolesol;
    
    if(status->params->boole_use_quadratic_model) {
        boolesol = tbkp_boolesol_quad_gurobi_get(status->instance, status->x, status->boole_grb_model);
    } else {
        boolesol = tbkp_boolesol_lin_gurobi_get(status->instance, status->x, status->boole_grb_model);
    }
    
    ++(status->stats->n_boole_called);
    status->stats->tot_time_boole += boolesol.time_to_compute;

    double local_lb = (double) (boolesol.lb_sum_profits + residual->sum_profits) *
                     (boolesol.lb_product_probabilities * residual->prod_probabilities);

    if(BB_VERBOSITY_CURRENT > BB_VERBOSITY_INFO) {
        printf("\t\tLB (Boole relaxation): %f (vs %f)\n", local_lb, status->solution->value);
    }

    if(current_node == 1u) {
        status->stats->boole_lb_at_root = local_lb;
        status->stats->time_to_compute_boole_at_root = boolesol.time_to_compute;
    }

    // Possibly update the best integer solution
    if(local_lb > status->stats->lb) {
        assert(local_lb > status->solution->value);

        status->stats->lb = local_lb;
        update_best_solution_from_boole(status, residual, &boolesol, local_lb);
    }

    tbkp_boolesol_free_inside(&boolesol);

    return local_lb;
}

/**
 * Solves the deterministic Knapsack Problem when all TB objects are fixed, at a leaf node.
 * If the new solution is better than the current best feasible solution, it updates it. 
 *
 * @param status                Current state of the algorithm.
 * @param residual              Current residual instance.
 */
static double solve_det_kp(TBKPBBAlgStatus* status, const TBKPBBResidualInstance *const residual) {
    size_t n_det_items = 0u;
    for(size_t i = 0; i < status->instance->n_items; ++i) {
        if(status->instance->probabilities[i] >= 1.0 - EPS) {
            // Deterministic items should all be unfixed.
            assert(status->x[i] == UNFIXED);
            n_det_items++;
        }
    }

    if(n_det_items > 0u) {
        cmb_item* cmb_items = malloc(n_det_items * sizeof(*cmb_items));
        cmb_stype sumW = 0;
        cmb_stype sumP = 0;

        if(!cmb_items) {
            printf("Cannot allocate memory for COMBO deterministic items\n");
            exit(EXIT_FAILURE);
        }

        size_t det_cnt = 0u;
        for(size_t i = 0; i < status->instance->n_items; ++i) {
            if((status->instance->probabilities[i] >= 1.0 - EPS)) {
                assert(status->x[i] == UNFIXED);

                cmb_items[det_cnt] = (cmb_item)
                {
                    .p = (cmb_itype) status->instance->profits[i],
                    .w = (cmb_itype) status->instance->weights[i],
                    .x = 0,
                    .pos = i
                };
                sumW += cmb_items[det_cnt].w;
                sumP += cmb_items[det_cnt].p;
                det_cnt++;
            }
        }

        cmb_stype myz = 0;

        // COMBO crashes if the total weight is less than the capacity!
        if((uint_fast32_t) sumW > residual->res_capacity) {
            myz = combo(&cmb_items[0], &cmb_items[n_det_items - 1], (cmb_stype)residual->res_capacity,
                    0, INT32_MAX, true, false);
        } else {
            myz = sumP;
            for(size_t i = 0u; i < n_det_items; ++i) {
                cmb_items[i].x = 1;
            }
        }

        // Compute the new solution's value
        uint_fast32_t sumprof = residual->sum_profits + (uint_fast32_t)myz;
        double zz = (double)sumprof * residual->prod_probabilities;

        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
                printf("\tCOMBO solution: %.3f\n", zz);
            }

        // Possibly update the incumbent
        if(zz > status->stats->lb) {
            assert(zz > status->solution->value);

            status->stats->lb = zz;
            status->solution->value = zz;
            status->solution->sum_profits = sumprof;
            status->solution->prod_probabilities = residual->prod_probabilities;

            if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
                printf("\tCOMBO found a new best feasible solution (%.3f vs %.3f)\n", zz, status->stats->lb);
            }

            // Time-bomb objects
            for(size_t i = 0u; i < status->instance->n_items; ++i) {
                if(status->x[i] != UNFIXED) {
                    status->solution->x[i] = status->x[i];
                }
            }

            // Non-time-bomb objects
            for(size_t i = 0u; i < n_det_items; ++i) {
                if(cmb_items[i].x) {
                    status->solution->x[cmb_items[i].pos] = true;
                }
            }
        }

        free(cmb_items); cmb_items = NULL;

        return zz;
    }

    return 0.0f;
}

/**
 * Branch on the current node and call `tbkp_bb_solve_node` on the two children nodes.
 *
 * @param status                Current algorithm state.
 * @param residual              Current residual instance.
 * @param parent_ub             UB at the current node (about to become father node).
 * @param current_node          Node number of the current node.
 * @param jbra                  Index of the item to branch on.
 */
static void branch(
        TBKPBBAlgStatus* status,
        TBKPBBResidualInstance residual,
        double parent_ub,
        size_t current_node,
        int jbra
) {
    // Return if the node limit has been reached
    if(status->params->max_nodes > 0u && *(status->n_nodes) >= status->params->max_nodes) {
        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("Maximum number of nodes (%zu) reached!\n", status->params->max_nodes);
        }

        // Use timeout to update the ub
        timeout(status, parent_ub, true);
        return;
    }

    // Left node: fix the item in the solution
    TBKPBBResidualInstance left_residual = {
            .prod_probabilities = residual.prod_probabilities * status->instance->probabilities[jbra],
            .sum_profits = residual.sum_profits + status->instance->profits[jbra],
            .res_capacity = residual.res_capacity - status->instance->weights[jbra]
    };

    // Fix the item: pack it!
    status->x[jbra] = FIXED_PACK;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Branching on item %d (profit: %" PRIuFAST32 ", weight %" PRIuFAST32 ", prob %.3f)\n",
                current_node, jbra, status->instance->profits[jbra], status->instance->weights[jbra],
                status->instance->probabilities[jbra]);
        printf("\tPacking the item\n");
        printf("\tResidual capacity in the child node: %" PRIuFAST32 "\n", left_residual.res_capacity);
        printf("\tSum of profits in the child node: %" PRIuFAST32 "\n", left_residual.sum_profits);
        printf("\tProduct of probabilities in the child node: %.3f\n", left_residual.prod_probabilities);
    }

    tbkp_bb_solve_node(status, parent_ub, left_residual, true /* early combo on 1-branch */);

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Returned from left node to father node\n", current_node);
    }

    // Return if the node limit has been reached
    if(status->params->max_nodes > 0u && *(status->n_nodes) >= status->params->max_nodes) {
        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("Maximum number of nodes (%zu) reached!\n", status->params->max_nodes);
        }

        // Use timeout to update the ub
        timeout(status, parent_ub, true);
        return;
    }

    // Right node: remove the item
    // The residual instance stays the same.

    // Fix the item: don't pack it!
    status->x[jbra] = FIXED_DONT_PACK;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Branching on item %d (profit: %" PRIuFAST32 ", weight %" PRIuFAST32 ", prob %.3f)\n",
               current_node, jbra, status->instance->profits[jbra], status->instance->weights[jbra],
               status->instance->probabilities[jbra]);
        printf("\tExcluding the item\n");
        printf("\tResidual capacity in the child node: %" PRIuFAST32 "\n", residual.res_capacity);
        printf("\tSum of profits in the child node: %" PRIuFAST32 "\n", residual.sum_profits);
        printf("\tProduct of probabilities in the child node: %.3f\n", residual.prod_probabilities);
    }

    tbkp_bb_solve_node(status, parent_ub, residual, false /* no early combo on 0-branch */);

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Returned from right node to father node\n", current_node);
    }

    // Unfix the item.
    status->x[jbra] = UNFIXED;
}

/** Solves a Branch-and-Bound node and possible updates the current solution.
 *
 * @param status      Current state of the algorithm.
 * @param parent_ub   Upper bound of the parent node.
 * @param residual    Residual instance at the current node.
 * @param early_combo Whether to get a LB from combo even before reaching a leaf node.
 */
static void tbkp_bb_solve_node(
        TBKPBBAlgStatus* status,
        double parent_ub,
        TBKPBBResidualInstance residual,
        _Bool early_combo)
{
    if(timeout(status, parent_ub, false)) {
        return;
    }

    (*(status->n_nodes))++;

    const size_t current_node = *(status->n_nodes);

    // Count the number of items that are still unfixed
    size_t n_unfixed_items = num_unfixed_items(status);

    // Define the items residual instance
    size_t* items = residual_instance(status, n_unfixed_items);

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Entering node\n", current_node);
        printf("\tCurrent sol = %.3f, residual capacity = %" PRIuFAST32 "\n",
                status->solution->value, residual.res_capacity);
        printf("\tSum of profits: %" PRIuFAST32 "\n", residual.sum_profits);
        printf("\tProduct of probabilities: %f\n", residual.prod_probabilities);
        printf("\tNumber of unfixed items: %zu\n", n_unfixed_items);
    }    

    double new_ub = DBL_MAX;
    double new_lb = 0.0;
    _Bool use_all_bounds = status->params->use_all_bounds_at_root && (current_node == 1u);

    /** START COMPUTATION OF NEW UBs AND/OR LBs AT THIS NODE. **/

    if(status->params->use_cr_bound || use_all_bounds) {
        CRBound cr_bound = get_cr_bound(status, current_node);

        if(cr_bound.local_ub < new_ub - EPS) {
            // New UB improved.
            new_ub = cr_bound.local_ub;
        }

        if(cr_bound.is_lb && cr_bound.local_ub > new_lb) {
            // New LB improved.
            new_lb = cr_bound.local_ub;
        }

        // Do not prune early if using all bounds, otherwise we cannot compute
        // the other bounds.
        if(cr_bound.should_prune && !use_all_bounds) {
            if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
                printf("[NODE %zu] Pruning node thanks to CR bound\n", current_node);
            }

            return;
        }
    }

    if(status->params->use_de_bounds || use_all_bounds) {
        DEBounds de_bounds = get_de_bounds(status, &residual, current_node, items, n_unfixed_items);

        if(de_bounds.local_ub < new_ub - EPS) {
            // Local UB improved.
            new_ub = de_bounds.local_ub;
        }

        if(de_bounds.local_lb > new_lb) {
            // Local LB improved.
            new_lb = de_bounds.local_lb;
        }

        // Do not prune early if using all bounds, otherwise we cannot compute
        // the other bounds.
        if(de_bounds.should_prune && !use_all_bounds) {
            if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
                printf("[NODE %zu] Pruning node thanks to DE bounds\n", current_node);
            }

            return;
        }
    }

    assert(current_node >= 1u);
    if((status->params->use_boole_bound && (current_node - 1u) % status->params->boole_bound_frequency == 0u) || use_all_bounds) {
        double boole_lb = get_boole_bound(status, &residual, current_node);

        if(boole_lb > new_lb) {
            // Local LB improved.
            new_lb = boole_lb;
        }
    }

    if((status->params->use_early_combo || use_all_bounds) && early_combo) {
        // It's smarter to solve the deterministic 01KP early, with the TB items fixed until now. When
        // branching on zero, we don't need to recompute this, because we inherit from the parent node.
        // We only need to call combo again on the 1 branches.

        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("[NODE %zu] Calling COMBO\n", current_node);
        }

        double cmb_lb = solve_det_kp(status, &residual);

        if(cmb_lb > new_lb) {
            // Local LB improved.
            new_lb = cmb_lb;
        }
    }

    /** END COMPUTATION OF NEW UBs AND/OR LBs AT THIS NODE. **/

    // Check if exploring this node gave a tighter UB than what we had.
    double local_ub = parent_ub;    
    if(local_ub == INITIAL_UB_PLACEHOLDER || new_ub < local_ub - EPS) {
        local_ub = new_ub;
    }

    // Nothing to check here. Just renaming new_lb to local_lb for consistency.
    const double local_lb = new_lb;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\tLocal LB for this node: %.6f\n", local_lb);
        printf("\tLocal UB for this node: %.6f\n", local_ub);
    }

    // Check if LB == UB; in this case, prune.
    if(local_lb >= local_ub - EPS && local_lb != INITIAL_LB_PLACEHOLDER && local_ub != INITIAL_UB_PLACEHOLDER) {
        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("[NODE %zu] LB and UB coincide (%.6f): closing node\n", current_node, local_lb);
        }

        free(items); items = NULL;
        return;
    }

    // Find the branching item
    int jbra = tbkp_bb_branch_item(status, &residual);

    if(jbra == NO_TIMEBOMB_ITEM_TO_BRANCH) {
        // All uncertain items have been fixed: solve a knapsack instance induced by the unfixed deterministic items.
        // Note: we are in a leaf node!

        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("[NODE %zu] No branching item: leaf node\n", current_node);
        }

        if(!status->params->use_early_combo) {
            // If not using early combo, we call combo only in the leaf nodes, i.e., here!
            solve_det_kp(status, &residual);
        }

        free(items); items = NULL;
        return;
    }

    branch(status, residual, local_ub, current_node, jbra);

    free(items); items = NULL;
}
/************************************************************
 * END OF LOCAL HELPER FUNCTIONS                            *
 ************************************************************/

TBKPBBSolution* tbkp_branch_and_bound(const TBKPInstance *const instance, TBKPBBStats* stats, TBKPParams* params) {
    TBKPBBSolution* solution = tbkp_sol_init(instance);
    TBKPBBFixedStatus* x = malloc(instance->n_items * sizeof(*x));

    if(!x) {
        printf("Cannot allocate memory for x variables in B&B!\n");
        exit(EXIT_FAILURE);
    }

    for(size_t t = 0u; t < instance->n_items; ++t) {
        x[t] = UNFIXED;
    }

    size_t n_nodes = 0u;

    GRBmodel* boole_grb_model = NULL;

    if(params->boole_use_quadratic_model) {
        boole_grb_model = tbkp_boolesol_quad_base_model(instance, params);
    } else {
        boole_grb_model = tbkp_boolesol_lin_base_model(instance, params);
    }

    TBKPBBAlgStatus status = {
        .instance = instance,
        .params = params,
        .solution = solution,
        .x = x,
        .stats = stats,
        .n_nodes = &n_nodes,
        .boole_grb_model = boole_grb_model
    };

    TBKPBBResidualInstance residual = {
        .prod_probabilities = 1.0,
        .sum_profits = 0u,
        .res_capacity = instance->capacity
    };

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("Starting solution at the root node\n");
    }

    tbkp_bb_solve_node(&status, INITIAL_UB_PLACEHOLDER, residual, true /* early combo */);

    free(x); x = NULL;

    stats->end_time = clock();
    stats->elapsed_time = (float)(stats->end_time - stats->start_time) / CLOCKS_PER_SEC;
    stats->n_nodes = n_nodes;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("Tree exploration completed after %zu nodes\n", n_nodes);
    }

    if(stats->ub == INITIAL_UB_PLACEHOLDER) {
        // UB was never updated or no open node left

        if(stats->elapsed_time < params->timeout && params->max_nodes == 0u) {
            // Instance solved to optimality
            stats->ub = stats->lb;
        } else {
            // Open instance
            stats->ub = FLT_MAX;
        }
    }

    stats->gap = (stats->ub - stats->lb) / stats->ub;

    GRBfreemodel(boole_grb_model);

    return solution;
}
