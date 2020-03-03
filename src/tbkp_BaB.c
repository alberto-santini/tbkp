//
// last update: Michele 03/03/2020.
//

#include "tbkp_BaB.h"
#include "tbkp_de_sol.h"
#include "utils/pdqsort_c.h"
#include <stdlib.h>
#include <stdio.h>
#include <combo.h>
#include <assert.h>

#ifndef BB_VERBOSITY_CURRENT
#define BB_VERBOSITY_CURRENT 1001
#endif

#define BB_VERBOSITY_INFO 1000
#define EPS 1e-6

TBKPSolution* tbkp_sol_init(const TBKPInstance *const instance) {
	TBKPSolution* solution = malloc(sizeof(*solution));

	if(!solution) {
		printf("Cannot allocate memory for solution\n");
		exit(EXIT_FAILURE);
	}

	solution->x = malloc(instance->n_items * sizeof(*solution->x));

	if(!solution->x) {
		printf("Cannot allocate memory for solution's x vector\n");
		exit(EXIT_FAILURE);
	}

	for(size_t i = 0u; i < instance->n_items; ++i) {
	    solution->x[i] = false;
	}

	solution->prod_probabilities = 1.0f;
	solution->sum_profits = 0u;
	solution->value = 0.0f;

	return solution;
}

void tbkp_sol_print(TBKPSolution* solution, const TBKPInstance *const instance) {
	printf("Solution with value %.2f (%" PRIuFAST32 " sum of profits, %.2f prod of probabilities)\n",
			solution->value, solution->sum_profits, solution->prod_probabilities);
	printf("Packed objects:\n");

	for(size_t i = 0u; i < instance->n_items; ++i) {
		if(solution->x[i]) {
			printf("\tObj %zu, profit %" PRIuFAST32 ", weight %"PRIuFAST32 ", prob %.2f\n",
					i, instance->profits[i], instance->weights[i], instance->probabilities[i]);
		}
	}
}

void tbkp_sol_free(TBKPSolution** solution_ptr) {
	free((*solution_ptr)->x); (*solution_ptr)->x = NULL;
	free(*solution_ptr); *solution_ptr = NULL;
}

TBKPSolution* tbkp_branch_and_bound(const TBKPInstance *const instance) {
	TBKPSolution* solution = tbkp_sol_init(instance);

	TBKPBBFixedStatus* x = malloc(instance->n_items * sizeof(*x));

	for(size_t t = 0u; t < instance->n_items; ++t) {
		x[t] = UNFIXED;
	}

	float prob_probabilities = 1.0f;
	uint_fast32_t sum_profits = 0;
	uint_fast32_t residual_capacity = instance->capacity;

	size_t nnodes = 0u;

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
		printf("[ROOT NODE] calling tbkp_bb_solve_node with capacity = %" PRIuFAST32 " and n_items = %" PRIuFAST32 "\n",
				instance->capacity, instance->n_items);

		for(size_t t = 0u; t < instance->n_items; ++t) {
			printf("[%zu]: profit %" PRIuFAST32 ", weight %" PRIuFAST32 ", probab %7.3lf\n",
					t, instance->profits[t], instance->weights[t], instance->probabilities[t]);
		}
		printf("\n");
	}

	tbkp_bb_solve_node(instance, &nnodes, x, prob_probabilities, sum_profits, residual_capacity, solution);

	free(x); x = NULL;
	return solution;
}

int tbkp_bb_branch_item(const TBKPInstance *const instance, const TBKPBBFixedStatus *const x, uint_fast32_t capacity) {
	// Find the first uncertain item that is not fixed and fits in the residual capacity
	for(size_t i = 0u; i < instance->n_items; ++i) {
		if(instance->probabilities[i] > 1.0 - EPS) continue;
		if(x[i] != UNFIXED) continue;
		if(instance->weights[i] > capacity) continue;

		return (int)i;
	}

	return NO_TIMEBOMB_ITEM_TO_BRANCH;
}


void tbkp_bb_solve_node(
        const TBKPInstance *const instance,
        size_t* nnodes,
        TBKPBBFixedStatus* x,
        float prod_probabilities,
        uint_fast32_t sum_profits,
        uint_fast32_t res_capacity,
        TBKPSolution *solution
) {
	(*nnodes)++;

    const size_t current_node = *nnodes;

	// Count the number of items that are still unfixed
	size_t n_unfixed_items = 0u;
	for(size_t j = 0u; j < instance->n_items; ++j) {
        if (x[j] == UNFIXED) {
            n_unfixed_items++;
        }
    }

	if(n_unfixed_items == 0) {
	    return;
	}

	// Define the residual instance
	size_t* items = malloc(n_unfixed_items * sizeof(*items));

	if(!items) {
	    printf("Cannot allocate memory for the residual instance's items\n");
	    exit(EXIT_FAILURE);
	}

	size_t cnt = 0;
	for(size_t j = 0; j < instance->n_items; j++) {
		if(x[j] < 0) {
			items[cnt++] = j;
		}
	}
	
	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("[NODE %zu] Best z = %.3f, residual capacity = %" PRIuFAST32 "\n",
	            current_node, solution->value, res_capacity);
		printf("\tSum of profits: %" PRIuFAST32 "\n", sum_profits);
		printf("\tProduct of probabilities: %f", prod_probabilities);
		printf("\tNumber of unfixed items: %zu", n_unfixed_items);
		printf("\tUnfixed items: ");
		for(size_t j = 0u; j < n_unfixed_items; ++j) {
		    printf("%zu ", items[j]);
		}
		printf("\n");
	}

	// Solve the deterministic relaxation
	TBKPDeterministicEqSol desol = tbkp_desol_get(instance, n_unfixed_items, items, res_capacity);
	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("Solution from the deterministic relaxation:\n");
	    tbkp_desol_print(&desol);
	}

	// Compute the local upper bound and possibly kill the node
	float localub = ((float)sum_profits + desol.ub) * prod_probabilities;
	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("UB from deterministic relaxation solution: %f\n", desol.ub);
	    printf("Sum of profits: %" PRIuFAST32 ", Product of probabilities: %f\n", sum_profits, prod_probabilities);
	    printf("Local UB: %f (vs %f)\n", localub, solution->value);
	}

	if(localub <= solution->value) {
		if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
		    printf("Node killed!\n");
		}

		free(items); items = NULL;
		return;
	}

	// Possibly update the incumbent
	if(desol.lb > solution->value ) {
		solution->value = desol.lb;

		for(size_t i = 0u; i < instance->n_items; ++i) {
			if(x[i] != UNFIXED) {
			    solution->x[i] = x[i];
			}
		}

		for(size_t i = 0u; i < desol.n_items; ++i) {
			solution->x[desol.items[i]] = true;
		}

		if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("Solution update: new value %f\n", solution->value);
        }
	}

	tbkp_desol_free_inside(&desol);

	// Find the branching item
	int jbra = tbkp_bb_branch_item(instance, x, res_capacity);

	if(jbra == NO_TIMEBOMB_ITEM_TO_BRANCH) {
		if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
		    printf("No branching item!\n");
		}

		// All uncertain items have been fixed: solve a knapsack instance induced by the unfixed deterministic items.
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
			    printf("Det combo solution: %ld, sumprof: %" PRIuFAST32 " => New z: %f (best z: %f)\n",
			            myz, sumprof, zz, solution->value);
			}

			// Possibly update the incumbent
			if(zz > solution->value ) {
				solution->value = zz;

				// Time-bomb objects
				for(size_t i = 0u; i < instance->n_items; ++i) {
					if(x[i] != UNFIXED) {
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
		}

		free(items); items = NULL;
		return;
	}

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

	tbkp_bb_solve_node(instance, nnodes, x, new_prod_probabilities, new_sum_profits, new_res_capacity, solution);

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("Returned from left node %zu to father node %zu\n", *nnodes, current_node);
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

	tbkp_bb_solve_node(instance, nnodes, x, new_prod_probabilities, new_sum_profits, new_res_capacity, solution);

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("Returned from right node %zu to father node %zu\n", *nnodes, current_node);
    }

	// Unfix the item.
	x[jbra] = UNFIXED;

	free(items); items = NULL;
}
