//
// last update: Michele 03/03/2020.
//

#include "tbkp_bb.h"
#include "tbkp_de_sol.h"
#include "utils/pdqsort_c.h"
#include "tbkp_boole_sol.h"
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
#define EPS 1e-6f

/************************************************************
 * LOCAL HELPER FUNCTIONS                                   *
 ************************************************************/

static void tbkp_bb_solve_node(const TBKPInstance* instance, size_t* nnodes, TBKPBBFixedStatus* x, float parent_ub,
        float prod_probabilities, uint_fast32_t sum_profits, uint_fast32_t res_capacity, _Bool early_combo,
        TBKPSolution *solution, TBKPBBStats* stats);

/** Finds the next time-bomb item to branch on. If no such item exists (because we either
 *  ran out of TB items, or there is no item which fits in the residual capacity), returns
 *  NO_TIMEBOMB_ITEM_TO_BRANCH.
 *
 * @param instance    The TBKP instance we are solving.
 * @param x           An array keeping track of the branching decisions:
 *                    x[i] == FIXED_DONT_PACK if item i has been fixed to NOT being packed
 *                    x[i] == FIXED_PACK if item i has been fixed to being packed
 *                    x[i] == UNFIXED if item i is unfixed.
 * @param capacity    Residual capacity of the knapsack.
 * @param sum_profits Sum-of-profits part of the objective function at current node.
 * @return            The index of the item to branch on, or NO_TIMEBOMB_ITEM_TO_BRANCH.
 */
static int tbkp_bb_branch_item(
        const TBKPInstance *const instance,
        const TBKPBBFixedStatus *const x,
        uint_fast32_t capacity,
        uint_fast32_t sum_profits)
{
    // Find the first uncertain item that is not fixed and fits in the residual capacity
    for(size_t i = 0u; i < instance->n_items; ++i) {
        if(instance->probabilities[i] > 1.0f - EPS) continue;
        if(x[i] != UNFIXED) continue;
        if(instance->weights[i] > capacity) continue;

        // Pruning: skip branch in case the solution value cannot increase
        float scorej = ((float)instance->profits[i] * instance->probabilities[i]) / (1.0f - instance->probabilities[i]);
        if(scorej < sum_profits) {
            continue;
        }

        return (int)i;
    }

    return NO_TIMEBOMB_ITEM_TO_BRANCH;
}

/**
 * Checks if timeout occurred while exploring the B&B tree.
 * It updates the TBKPBBStats object with
 * @param stats     Solution statistics
 * @param parent_ub UB of the parent node of the current node.
 * @param nnodes    Number of nodes being solved.
 * @return          Returns true iff a timeout occurred.
 */
static _Bool timeout(
        TBKPBBStats *const stats,
        float parent_ub,
        size_t nnodes
) {
	clock_t current_time = clock();
    float el_time = (float)(current_time - stats->start_time) / CLOCKS_PER_SEC;

    if(el_time > stats->timeout) {
        stats->n_nodes = nnodes;

        // Upon timeout, the UB is the worst UB of the open nodes.
        if(stats->ub == INITIAL_UB_PLACEHOLDER || (parent_ub != INITIAL_UB_PLACEHOLDER && parent_ub > stats->ub)) {
            stats->ub = parent_ub;
        }

        return true;
    }

	return false;
}

/**
 * Returns the number of unfixed TB items at current node.
 * @param instance  TBKP instance.
 * @param x         Array with the status of each item.
 * @return          The number of items with UNFIXED status.
 */
static size_t num_unfixed_items(
        const TBKPInstance *const instance,
        const TBKPBBFixedStatus *const x
) {
	size_t n_unfixed_items = 0u;
	for(size_t j = 0u; j < instance->n_items; ++j) {
        if (x[j] == UNFIXED) {
            n_unfixed_items++;
        }
    }
	return n_unfixed_items;
}

/**
 * Gives the residual instance at the current node, i.e., the instance without
 * all fixed TB items.
 *
 * @param instance          TBKP instance.
 * @param x                 List of item statuses.
 * @param n_unfixed_items   Number of unfixed items in list x.
 * @return                  An array with the indices of items in the residual instance.
 */
static size_t* residual_instance(
        const TBKPInstance *const instance,
        const TBKPBBFixedStatus *const x,
        size_t n_unfixed_items)
{
	size_t* items = malloc(n_unfixed_items * sizeof(*items));
	if(!items) {
	    printf("Cannot allocate memory for the residual instance's items\n");
	    exit(EXIT_FAILURE);
	}

	size_t cnt = 0;
	for(size_t j = 0; j < instance->n_items; j++) {
		if(x[j] == UNFIXED) {
			items[cnt++] = j;
		}
	}

	return items;
}

/**
 * Structure to contain info from the Deterministic Equivalent bounds,
 * which are useful in the B&B nodes. It contains the UB and LB at the
 * local(current) node, plus a flag indicating if the current node should
 * be prunned - i.e., because the local UB is worse than the best-known
 * LB.
 */
typedef struct {
	float local_ub;
	float local_lb;
	_Bool should_prune;
} DEBounds;

/**
 * Updates the current best solution from the LB-solution gotten from the
 * Deterministic Equivalent solution.
 *
 * @param instance              TBKP instance.
 * @param solution              Current solution (to update).
 * @param desol                 Deterministic Equivalent solution.
 * @param x                     List of items statuses.
 * @param local_lb              Local LB, i.e., obj value of the DE solution.
 * @param prod_probabilities    Product-of-probabilites component of the obj function.
 * @param sum_profits           Sum-of-profits component of the obj function.
 */
static void update_best_solution_from_de(
        const TBKPInstance *const instance,
        TBKPSolution* solution,
        const TBKPDeterministicEqSol *const desol,
        const TBKPBBFixedStatus *const x,
        float local_lb,
        float prod_probabilities,
        uint_fast32_t  sum_profits
) {
    solution->value = local_lb;
    solution->prod_probabilities = desol->lb_product_probabilities * prod_probabilities;
    solution->sum_profits = desol->lb_sum_profits + sum_profits;

    for(size_t i = 0u; i < instance->n_items; ++i) {
        if(x[i] != UNFIXED) {
            solution->x[i] = x[i];
        }
    }

    for(size_t i = 0u; i < desol->n_items; ++i) {
        solution->x[desol->items[i]] = true;
    }

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("\tSolution update from DE: new value %f\n", solution->value);
    }
}

/**
 * Updates the current best solution from the LB-solution gotten from the
 * Boole bound solution.
 *
 * @param instance              TBKP instance.
 * @param solution              Current solution (to update).
 * @param desol                 Boole bound solution.
 * @param x                     List of items statuses.
 * @param local_lb              Local LB, i.e., obj value of the Boole solution.
 * @param prod_probabilities    Product-of-probabilites component of the obj function.
 * @param sum_profits           Sum-of-profits component of the obj function.
 */
static void update_best_solution_from_boole(
        const TBKPInstance *const instance,
        TBKPSolution* solution,
        const TBKPBooleSol *const boolesol,
        const TBKPBBFixedStatus *const x,
        float local_lb,
        float prod_probabilities,
        uint_fast32_t  sum_profits
) {
    solution->value = local_lb;
    solution->prod_probabilities = boolesol->lb_product_probabilities * prod_probabilities;
    solution->sum_profits = boolesol->lb_sum_profits + sum_profits;

    for(size_t i = 0u; i < instance->n_items; ++i) {
        if(x[i] != UNFIXED) {
            solution->x[i] = x[i];
        }
    }

    for(size_t i = 0u; i < boolesol->n_items; ++i) {
        solution->x[boolesol->items[i]] = true;
    }

    if(BB_VERBOSITY_CURRENT > BB_VERBOSITY_INFO) {
        printf("\tSolution update from Boole: new value %f\n", solution->value);
    }
}

/**
 * Gets the Deterministic Equivalent bounds (UB and LB) at the current node.
 * If the UB is worse than the current best LB, it signal to prune the current node.
 * It the LB is better than the current best LB, it updates it.
 *
 * @param instance              TBKP instance.
 * @param stats                 Solution stats.
 * @param solution              Current solution.
 * @param x                     List with item statuses.
 * @param items                 List with the indices of items in the current residual instance.
 * @param n_unfixed_items       Number of unfixed items in the current residual instance.
 * @param prod_probabilities    Product-of-probabilities component of the objective function.
 * @param sum_profits           Sum-of-profits component of the objective function.
 * @param res_capacity          Residual knapsack capacity at the current node.
 * @return                      The DEBounds object containing the DE bounds and the pruning flag.
 */
static DEBounds get_de_bounds(
        size_t current_node,
        const TBKPInstance *const instance,
        TBKPBBStats* stats,
        TBKPSolution* solution,
        const TBKPBBFixedStatus *const x,
        size_t* items,
        size_t n_unfixed_items,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity
) {
	// Solve the deterministic relaxation
	TBKPDeterministicEqSol desol = tbkp_desol_get(instance, n_unfixed_items, items, res_capacity);

	// Compute the local upper bound and possibly kill the node
	float local_ub = ((float)sum_profits + desol.ub) * prod_probabilities;
	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("\tUB from deterministic relaxation solution: %f\n", desol.ub);
	    printf("\tSum of profits: %" PRIuFAST32 ", Product of probabilities: %f\n", sum_profits, prod_probabilities);
	    printf("\tLocal UB: %f (vs %f)\n", local_ub, solution->value);
	}

	if(local_ub <= solution->value) {
        tbkp_desol_free_inside(&desol);
		free(items); items = NULL;
		return (DEBounds){.local_ub = local_ub, .local_lb = 0.0f, .should_prune = true};
	}

	float local_lb = (float) (desol.lb_sum_profits + sum_profits) *
                     (desol.lb_product_probabilities * prod_probabilities);

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("\tLB from deterministic relaxation solution: %f\n", desol.lb);
	    printf("\tSum of profits: %d (fixed %" PRIuFAST32 "), ", (int) desol.lb_sum_profits, sum_profits);
	    printf("\tProduct of probabilities %f (fixed %f)\n", desol.lb_product_probabilities, prod_probabilities);
	    printf("\tLocal LB: %f (vs %f)\n", local_lb, solution->value);
	}

    if(current_node == 1u) {
        stats->de_lb_at_root = local_lb;
        stats->de_ub_at_root = local_ub;
        stats->time_to_compute_de_at_root = desol.time_to_compute;
    }

	if(local_lb > stats->lb) {
	    stats->lb = local_lb;
	}

    // Possibly update the incumbent
	if(local_lb > solution->value ) {
        update_best_solution_from_de(instance, solution, &desol, x, local_lb, prod_probabilities, sum_profits);
	}

	tbkp_desol_free_inside(&desol);

	return (DEBounds){.local_ub = local_ub, .local_lb = local_lb, .should_prune = false};
}

/**
 * Gets the Boole bound (LB) at the current node.
 * It the LB is better than the current best LB, it updates it.
 *
 * @param instance              TBKP instance.
 * @param stats                 Solution stats.
 * @param solution              Current solution.
 * @param x                     List with item statuses.
 * @param items                 List with the indices of items in the current residual instance.
 * @param n_unfixed_items       Number of unfixed items in the current residual instance.
 * @param prod_probabilities    Product-of-probabilities component of the objective function.
 * @param sum_profits           Sum-of-profits component of the objective function.
 * @param res_capacity          Residual knapsack capacity at the current node.
 * @return                      The DEBounds object containing the DE bounds and the pruning flag.
 */
float get_boole_bound(
        size_t current_node,
        const TBKPInstance *const instance,
        TBKPBBStats* stats,
        TBKPSolution* solution,
        const TBKPBBFixedStatus *const x,
        size_t* items,
        size_t n_unfixed_items,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity
) {
    TBKPBooleSol boolesol = tbkp_boolesol_lin_gurobi_get(instance , n_unfixed_items, items, res_capacity);

    float local_lb = (float) (boolesol.lb_sum_profits + sum_profits) *
                     (boolesol.lb_product_probabilities * prod_probabilities);

    if(BB_VERBOSITY_CURRENT > BB_VERBOSITY_INFO) {
        printf("\tLB from Boole relaxation solution: %.3f\n", boolesol.lb);
        printf("\tSum of profits: %d (fixed %" PRIuFAST32 ")\n", (int) boolesol.lb_sum_profits, sum_profits);
        printf("\tProduct of probabilities %f (fixed %f)\n", boolesol.lb_product_probabilities, prod_probabilities);
        printf("\tLocal LB: %f (vs %f)\n", local_lb, solution->value);
    }

    if(current_node == 1u) {
        stats->boole_lb_at_root = local_lb;
        stats->time_to_compute_boole_at_root = boolesol.time_to_compute;
    }

    if(local_lb > stats->lb) {
        stats->lb = local_lb;
    }

    // Possibly update the incumbent
    if(local_lb > solution->value) {
        update_best_solution_from_boole(instance, solution, &boolesol, x, local_lb, prod_probabilities, sum_profits);
    }

    tbkp_boolesol_free_inside(&boolesol);

    return local_lb;
}

/**
 * Solves the deterministic Knapsack Problem when all TB objects are fixed, at a leaf node.
 * If the new solution is better than `solution`, it updates it.
 *
 * @param instance              TBKP instance.
 * @param stats                 Solution statistics.
 * @param solution              Current solution.
 * @param x                     Array of item statuses.
 * @param prod_probabilities    Product-of-probabilities component of the obj function.
 * @param sum_profits           Sum-of-profits component of the obj function.
 * @param res_capacity          Residual knapsack capacity of the current instance.
 */
static void solve_det_kp(
        const TBKPInstance *const instance,
        TBKPBBStats* stats,
        TBKPSolution* solution,
        const TBKPBBFixedStatus *const x,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity
) {
		size_t n_det_items = 0u;
		for(size_t i = 0; i < instance->n_items; ++i) {
		    if(instance->probabilities[i] > 1.0 - EPS) {
		        // Deterministic items should all be unfixed.
		        assert(x[i] == UNFIXED);
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
			for(size_t i = 0; i < instance->n_items; ++i) {
				if((instance->probabilities[i] > 1.0 - EPS)) {
				    assert(x[i] == UNFIXED);

					cmb_items[det_cnt] = (cmb_item)
					{
						.p = (cmb_itype) instance->profits[i],
						.w = (cmb_itype) instance->weights[i],
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
			if((uint_fast32_t) sumW > res_capacity) {
				myz = combo(&cmb_items[0], &cmb_items[n_det_items - 1], (cmb_stype)res_capacity, 0, INT32_MAX, true, false);
			} else {
				myz = sumP;
				for(size_t i = 0u; i < n_det_items; ++i) {
				    cmb_items[i].x = 1;
				}
			}

			// Compute the new solution value and possibly update the incumbent
			uint_fast32_t sumprof = sum_profits + (uint_fast32_t)myz;
			float zz = (float)sumprof * prod_probabilities;

			if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
			    printf("\tCOMBO solution: %ld, sumprof: %" PRIuFAST32 " => New z: %f (best z: %f)\n",
			            myz, sumprof, zz, solution->value);
			}

			if(zz > stats->lb) {
			    stats->lb = zz;
			}

			// Possibly update the incumbent
			if(zz > solution->value) {
				solution->value = zz;

				// Time-bomb objects
				for(size_t i = 0u; i < instance->n_items; ++i) {
					if(x[i] == FIXED_PACK) {
                        solution->x[i] = x[i];
                    }
				}

				// Non-time-bomb objects
				for(size_t i = 0u; i < n_det_items; ++i) {
					if(cmb_items[i].x) {
						solution->x[cmb_items[i].pos] = true;
					}
				}
			}

            free(cmb_items); cmb_items = NULL;
		}
}

/**
 * Branch on the current node and call `tbkp_bb_solve_node` on the two children nodes.
 *
 * @param instance              TBKP instance.
 * @param nnodes                Number of nodes counter.
 * @param x                     Array with item statuses.
 * @param local_ub              UB at the current node.
 * @param prod_probabilities    Product-of-probabilities component of the obj function.
 * @param sum_profits           Sum-of-profits component of the obj function.
 * @param res_capacity          Residual knapsack capacity at the current node.
 * @param solution              Current solution.
 * @param stats                 Solution stats.
 * @param current_node          Node number of the current node.
 * @param jbra                  Index of the item to branch on.
 */
static void branch(
        const TBKPInstance *const instance,
        size_t* nnodes,
        TBKPBBFixedStatus* x,
        float local_ub,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity,
        TBKPSolution* solution,
        TBKPBBStats* stats,
        size_t current_node,
        int jbra
) {
	// Left node: fix the item in the solution
	uint_fast32_t new_res_capacity = res_capacity - instance->weights[jbra];
	uint_fast32_t new_sum_profits = sum_profits + instance->profits[jbra];
	float new_prod_probabilities = prod_probabilities * instance->probabilities[jbra];

	// Fix the item: pack it!
	x[jbra] = FIXED_PACK;

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("[NODE %zu] Branching on item %d (profit: %" PRIuFAST32 ", weight %" PRIuFAST32 ", prob %.3f)\n",
	            current_node, jbra, instance->profits[jbra], instance->weights[jbra], instance->probabilities[jbra]);
	    printf("\tPacking the item\n");
	    printf("\tResidual capacity in the child node: %" PRIuFAST32 "\n", new_res_capacity);
	    printf("\tSum of profits in the child node: %" PRIuFAST32 "\n", new_sum_profits);
	    printf("\tProduct of probabilities in the child node: %.3f\n", new_prod_probabilities);
	}

	tbkp_bb_solve_node(
	        instance, nnodes, x, local_ub, new_prod_probabilities, new_sum_profits, new_res_capacity,
	        true /* early combo on 1-branch */, solution, stats);

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("[NODE %zu] Returned from left node to father node\n", current_node);
	}

	// Right node: remove the item
	new_res_capacity = res_capacity;
    new_sum_profits = sum_profits;
    new_prod_probabilities = prod_probabilities;

    // Fix the item: don't pack it!
	x[jbra] = FIXED_DONT_PACK;

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Branching on item %d (profit: %" PRIuFAST32 ", weight %" PRIuFAST32 ", prob %.3f)\n",
               current_node, jbra, instance->profits[jbra], instance->weights[jbra], instance->probabilities[jbra]);
        printf("\tExcluding the item\n");
        printf("\tResidual capacity in the child node: %" PRIuFAST32 "\n", new_res_capacity);
        printf("\tSum of profits in the child node: %" PRIuFAST32 "\n", new_sum_profits);
        printf("\tProduct of probabilities in the child node: %.3f\n", new_prod_probabilities);
	}

	tbkp_bb_solve_node(
	        instance, nnodes, x, local_ub, new_prod_probabilities, new_sum_profits, new_res_capacity,
	        false /* no early combo on 0-branch */, solution, stats);

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("[NODE %zu] Returned from right node to father node\n", current_node);
    }

	// Unfix the item.
	x[jbra] = UNFIXED;
}

/** Solves a Branch-and-Bound node and possible updates the current solution.
 *
 * @param instance              TBKP instance we are solving.
 * @param nnodes                Number of nodes in the BB tree.
 * @param x                     Vector keeping track of fixed/unfixed tb items.
 * @param parent_ub             Upper Bound of the parent node.
 * @param prod_probabilities    Product-of-probabilities part of the objective function.
 * @param sum_profits           Sum-of-profits part of the objective function.
 * @param res_capacity          Residual capacity of the knapsack at this node.
 * @param early_combo           If true, solve a deterministic knapsack with the TB items
 *                              currently fixed, to obtain a lower bound.
 * @param solution              Current solution.
 * @param stats                 Solution statistics.
 */
static void tbkp_bb_solve_node(
        const TBKPInstance *const instance,
        size_t* nnodes,
        TBKPBBFixedStatus* x,
        float parent_ub,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity,
        _Bool early_combo,
        TBKPSolution *solution,
        TBKPBBStats* stats
) {
    if(timeout(stats, parent_ub, *nnodes)) {
		return;
	}

	(*nnodes)++;

    const size_t current_node = *nnodes;

	// Count the number of items that are still unfixed
	size_t n_unfixed_items = num_unfixed_items(instance, x);
	if(n_unfixed_items == 0) {
	    return;
	}

	// Define the residual instance
	size_t* items = residual_instance(instance, x, n_unfixed_items);

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("[NODE %zu] Entering node\n", current_node);
	    printf("\tCurrent sol = %.3f, residual capacity = %" PRIuFAST32 "\n", solution->value, res_capacity);
		printf("\tSum of profits: %" PRIuFAST32 "\n", sum_profits);
		printf("\tProduct of probabilities: %f\n", prod_probabilities);
		printf("\tNumber of unfixed items: %zu\n", n_unfixed_items);
	}	

	float local_ub = parent_ub;
	float local_lb = 0.0f;

	if(stats->use_de_bounds) {
		DEBounds de_bounds = get_de_bounds(current_node,
		        instance, stats, solution, x, items, n_unfixed_items, prod_probabilities, sum_profits, res_capacity);

        if(local_ub == INITIAL_UB_PLACEHOLDER || de_bounds.local_ub < local_ub) { local_ub = de_bounds.local_ub; }
        if(de_bounds.local_lb > local_lb) { local_lb = de_bounds.local_lb; }

		if(de_bounds.should_prune) {
			if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
				printf("[NODE %zu] Pruning node thanks to DE bounds.\n", current_node);
			}

			return;
		}
	}

	if(stats->use_boole_bound && current_node % stats->boole_bound_frequency == 1u) {
		float boole_lb = get_boole_bound(current_node,
		        instance, stats, solution, x, items, n_unfixed_items, prod_probabilities, sum_profits, res_capacity);

		if(boole_lb > local_lb) { local_lb = boole_lb; }
	}

	if(local_lb > local_ub - EPS && local_lb > EPS && local_ub > EPS) {
	    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	        printf("[NODE %zu] LB and UB coincide (%.3f). Closing node.\n", current_node, local_lb);
	    }

	    free(items); items = NULL;
	    return;
	}

	if(stats->use_early_combo && early_combo) {
		// It's smarter to solve the deterministic 01KP early, with the TB items fixed until now. When
		// branching on zero, we don't need to recompute this, because we inherit from the parent not.
		// We only need to call combo again on the 1 branches.

		if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
		    printf("[NODE %zu] Calling COMBO.\n", current_node);
		}

		solve_det_kp(instance, stats, solution, x, prod_probabilities, sum_profits, res_capacity);
	}

    // Find the branching item
    int jbra = tbkp_bb_branch_item(instance, x, res_capacity, sum_profits);

    if(jbra == NO_TIMEBOMB_ITEM_TO_BRANCH) {
        // All uncertain items have been fixed: solve a knapsack instance induced by the unfixed deterministic items.
        // Note: we are in a leaf node!

        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("[NODE %zu] No branching item. Leaf node.\n", current_node);
        }

        if(!stats->use_early_combo) {
            // If not using early combo, we call combo only in the leaf nodes,
            // i.e., here!
            solve_det_kp(instance, stats, solution, x, prod_probabilities, sum_profits, res_capacity);
        }

        free(items); items = NULL;
        return;
    }

	branch(
	        instance, nnodes, x, local_ub, prod_probabilities, sum_profits, res_capacity,
	        solution, stats, current_node, jbra);

	free(items); items = NULL;
}
/************************************************************
 * END OF LOCAL HELPER FUNCTIONS                            *
 ************************************************************/

TBKPSolution* tbkp_branch_and_bound(const TBKPInstance *const instance, TBKPBBStats* stats) {
    TBKPSolution* solution = tbkp_sol_init(instance);
    TBKPBBFixedStatus* x = malloc(instance->n_items * sizeof(*x));

    if(!x) {
        printf("Cannot allocate memory for x variables in B&B!\n");
        exit(EXIT_FAILURE);
    }

    for(size_t t = 0u; t < instance->n_items; ++t) {
        x[t] = UNFIXED;
    }

    float prod_probabilities = 1.0f;
    uint_fast32_t sum_profits = 0;
    uint_fast32_t residual_capacity = instance->capacity;

    size_t nnodes = 0u;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("Starting solution at the root node.\n");
    }

    tbkp_bb_solve_node(instance, &nnodes, x, INITIAL_UB_PLACEHOLDER, prod_probabilities, sum_profits,
                       residual_capacity, true /* early combo */, solution, stats);

    free(x); x = NULL;

    stats->end_time = clock();
    stats->elapsed_time = (float)(stats->end_time - stats->start_time) / CLOCKS_PER_SEC;
    stats->n_nodes = nnodes;

    if(stats->ub == INITIAL_UB_PLACEHOLDER) {
        // UB was never updated

        if(stats->elapsed_time < stats->timeout) {
            // Instance solved to optimality
            stats->ub = stats->lb;
        }
    }

    stats->gap = (stats->ub - stats->lb) / stats->ub;

    return solution;
}