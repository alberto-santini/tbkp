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


/************michele: rilassamento continuo (inizio) ***************/
#define EPSC 1e-10

/**********************************************************************/
double bin_search(int n, double *x, double *d, uint_fast32_t *p, double *q)
/**********************************************************************/
{
// next point is z = x + lambda *x   with lambda \in [0,1]
        double bestobj = 0.0;
        double bestl = -1.0;

        double ptot = 0.0;
        double product = 1.0;

        // iterative loop
        double l0 = 0.0;
        double l1 = 1.0;
        while ( l1 - l0 > EPSC )
        {
                // define two lambda values in the current interval
                double lgap = l1 - l0;
                double delta = lgap / 3;
                double beta = l0 + delta;
                double gamma = l0 + 2 * delta;

                // evaluate the first point
                ptot = 0.0;
                product = 1.0;
                for ( int j = 0; j < n; j++ )
                {
                        double zj = x[j] + beta * d[j];
                        ptot += p[j] * zj;
                        double alpha = 1.0 - q[j] * zj;
                        product = product * alpha;
                }
                double objbeta = ptot * product;
                if ( objbeta > bestobj )
                {
                        bestobj = objbeta;
                        bestl = beta;
                }

                // evaluate the second point
                ptot = 0.0;
                product = 1.0;
                for ( int j = 0; j < n; j++ )
                {
                        double zj = x[j] + gamma * d[j];
                        ptot += p[j] * zj;
                        double alpha = 1.0 - q[j] * zj;
                        product = product * alpha;
                }
                double objgamma = ptot * product;
                if ( objgamma > bestobj )
                {
                        bestobj = objgamma;
                        bestl = gamma;
                }

                if ( objbeta < objgamma )
                {
                        // remove interval [l0, beta]
                        l0 = beta;
                }
                else
                {
                        // remove interval [gamma, l1]
                        l1 = gamma;
                }
        }

	return bestl;
}


/**********************************************************************/
void solve_Dantzig(int n, double *p, uint_fast32_t *w, int c, TBKPBBFixedStatus* xbra, double *y)
/**********************************************************************/
{
        uint_fast32_t cres = c;
        int *flag = (int *) calloc(n, sizeof(int));
        for ( int j = 0; j < n; j++ )
        {
                y[j] = 0;
                flag[j] = 0;
		// handle items that have been fixed by branching
		if ( xbra[j] != -1 ) 
		{
			flag[j] = 1;
			if ( xbra[j] == 1 ) 
			{
				y[j] = 1.0;
				cres -= w[j];
			}
		}
        }
        for ( int cont = 0; cont < n; cont++ )
        {
                if ( cres <= 0 ) break;

                double maxratio = 0.0;
                int index = -1;
                for ( int j = 0; j < n; j++ )
                {
                        if ( flag[j] == 0 )
                        {
                                double ratio = p[j]/w[j];
                                if ( ratio > maxratio )
                                {
                                        maxratio = ratio;
                                        index = j;
                                }
                        }
                }

                if ( index < 0 ) break;
                flag[index] = 1;
                if ( cres > w[index] )
                {
                        y[index] = 1.0;
                        cres -= w[index];
                }
                else
                {
                        y[index] = 1.0 * cres / w[index];
                        cres = 0;
                }
        }
        free(flag);

}



/**********************************************************************/
float solve_cont(const TBKPInstance *const instance, TBKPBBFixedStatus* xbra, double *xc)
/**********************************************************************/
{
	int n = instance->n_items;
	double *a = (double *) calloc(n, sizeof(double)); // \alpha variables
	double *q = (double *) calloc(n, sizeof(double)); // q probabilities
	for ( int j = 0; j < n; j++ ) 
	{
		if ( xbra[j] == 1 ) xc[j] = 1.0; else xc[j] = 0.0;
		a[j] = 1.0;
		q[j] = 1.0 - instance->probabilities[j];
	}
	double product = 1.0;
	double ptot = 0.0;
	double *g = (double *) calloc(n, sizeof(double)); // gradient
	double *d = (double *) calloc(n, sizeof(double)); // improving direction (if any)
	double *y = (double *) calloc(n, sizeof(double)); // next solution
	float objval = 0.0;
	while ( 1 ) 
	{
		float oldobj = objval;

		// compute the gradient
		product = 1.0;
                for ( int j = 0; j < n; j++ ) product *= (1.0 - q[j]*xc[j]);
                for ( int j = 0; j < n; j++ )
                {
                        double value = instance->profits[j] - ptot*q[j]/(1.0 - q[j]*xc[j]);
                        g[j] = product * value;
                }
		
		// find a candidate direction for improving
		solve_Dantzig(n, g, instance->weights, instance->capacity, xbra, y);
		for ( int j = 0; j < n; j++ ) d[j] = y[j] - xc[j];

		// check if the direction is improving
		double delta = 0.0;
		for ( int j = 0; j < n; j++ ) delta += g[j] * d[j];
		if ( delta < EPSC ) break;

		// determine the optimal value for parameter \lambda
		double lambda = bin_search(n, xc, d, instance->profits, q);

		// determine the next point
		ptot = 0.0;
		product = 1.0;
		for ( int j = 0; j < n; j++ )
                {
                        double zj = xc[j] +  lambda * d[j];
                        ptot += instance->profits[j] * zj;
                        a[j] = 1.0 - q[j] * zj;
                        product = product * a[j];
                        xc[j] = zj;
                }
                objval = ptot * product;
		
		// break if no improvement
		if ( objval - oldobj < EPSC ) break;
	}

	free(y);
	free(d);
	free(g);
	free(q);
	free(a);

	return objval;
}

/************michele: rilassamento continuo (fine) ***************/


/************************************************************
 * LOCAL HELPER FUNCTIONS                                   *
 ************************************************************/

static void tbkp_bb_solve_node(
        TBKPBBAlgStatus* status, float parent_ub, TBKPBBResidualInstance residual, _Bool early_combo, size_t max_nodes);

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
        if(status->instance->probabilities[i] > 1.0f - EPS) continue;
        if(status->x[i] != UNFIXED) continue;
        if(status->instance->weights[i] > residual->res_capacity) continue;

        // Pruning: skip branch in case the solution value cannot increase
        float scorej = ((float)status->instance->profits[i] * status->instance->probabilities[i]) /
                (1.0f - status->instance->probabilities[i]);

        if(scorej < residual->sum_profits) {
            continue;
        }

        return (int)i;
    }

    return NO_TIMEBOMB_ITEM_TO_BRANCH;
}

/**
 * Checks if timeout occurred while exploring the B&B tree.
 * It updates the TBKPBBStats object with
 * @param status    Current state of the algorithm.
 * @param parent_ub Upper bound at the parent node.
 * @return          Returns true iff a timeout occurred.
 */
static _Bool timeout(TBKPBBAlgStatus* status, float parent_ub) {
	clock_t current_time = clock();
    float el_time = (float)(current_time - status->stats->start_time) / CLOCKS_PER_SEC;

    if(el_time > status->params->timeout) {
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
 * local(current) node, plus a flag indicating if the current node should
 * be prunned - i.e., because the local UB is worse than the best-known
 * LB.
 */
typedef struct {
	float local_ub;
	float local_lb;
	_Bool should_prune;
	float local_ub2;	// ub from continuous relaxation
} DEBounds;

/**
 * Updates the current best solution from the LB-solution gotten from the
 * Deterministic Equivalent solution.
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
        float local_lb
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
        printf("\tSolution update from DE: new value %f\n", status->solution->value);
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
        float local_lb
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
        printf("\tSolution update from Boole: new value %f\n", status->solution->value);
    }
}

/**
 * Gets the Deterministic Equivalent bounds (UB and LB) at the current node.
 * If the UB is worse than the current best LB, it signal to prune the current node.
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

	// Compute the local upper bound and possibly kill the node
	float local_ub = ((float)residual->sum_profits + desol.ub) * residual->prod_probabilities;
	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("\tUB from deterministic relaxation solution: %f\n", desol.ub);
	    printf("\tSum of profits: %" PRIuFAST32 ", Product of probabilities: %f\n",
	            residual->sum_profits, residual->prod_probabilities);
	    printf("\tLocal UB: %f (vs %f)\n", local_ub, status->solution->value);
	}

	if(local_ub <= status->solution->value) {
        tbkp_desol_free_inside(&desol);
		free(items); items = NULL;
		return (DEBounds){.local_ub = local_ub, .local_lb = 0.0f, .should_prune = true};
	}

	// compute the continuous relaxation upper bound
	double *xc = (double *) calloc(status->instance->n_items, sizeof(double));
	float local_ub2 = solve_cont(status->instance, status->x, xc);
	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("\tUB from continuous relaxation: %f (vs %f)\n", local_ub2, status->solution->value);
	}
	if(local_ub2 <= status->solution->value) {
		free(items); items = NULL;
		return (DEBounds){.local_ub2 = local_ub2, .local_lb = 0.0f, .should_prune = true};
	}
	// check if the solution of the continuous relaxation is integer
	int integer_sol = 1;
	for(size_t i = 0u; i < status->instance->n_items; ++i) {	
		if ( (xc[i] > EPS) && (xc[i] < 1.0-EPS) )
		{
			integer_sol = 0;
			break;
		}
	}
///////////////////////Alberto
	if ( integer_sol ) 	
	{
		float lb_val = local_ub2;
		if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
			printf("\n*** WARNING: continuous relaxation has an integer solution with value %f ***\n", lb_val);
			printf("*** possibly update the incumbent and execute a backtracking ***\n\n");
		}
	}
///////////////////////Alberto

	
	float local_lb = (float) (desol.lb_sum_profits + residual->sum_profits) *
                     (desol.lb_product_probabilities * residual->prod_probabilities);


///////////////////////Alberto
	if ( desol.lb > local_lb ) 
	{
		printf("\n\n*** WARNING: current lower bound is %f (instead of %f) ***\n\n", desol.lb, local_lb);
	}
///////////////////////Alberto


	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("\tLB from deterministic relaxation solution: %f\n", desol.lb);
	    printf("\tSum of profits: %" PRIuFAST32 " (fixed %" PRIuFAST32 "), ",
	            desol.lb_sum_profits, residual->sum_profits);
	    printf("Product of probabilities %f (fixed %f)\n",
	            desol.lb_product_probabilities, residual->prod_probabilities);
	    printf("\tLocal LB: %f (vs %f)\n", local_lb, status->solution->value);
	}


///////////////////////Alberto
	if ( (desol.lb >= local_ub) || (desol.lb >= local_ub2) || (local_lb >= local_ub) || (local_lb >= local_ub2) ) 
	{
		printf("\n\n*** WARNING: local lower bounds are %f and %f ***\n", desol.lb, local_lb);
		printf("*** WARNING: local upper bounds are %f and %f ***\n", local_ub, local_ub2);
		printf("*** WARNING: kill the node ***\n\n");
	}
///////////////////////Alberto



    if(current_node == 1u) {
        status->stats->de_lb_at_root = local_lb;
        status->stats->de_ub_at_root = local_ub;
        status->stats->time_to_compute_de_at_root = desol.time_to_compute;
    }

	if(local_lb > status->stats->lb) {
	    status->stats->lb = local_lb;
	}

    // Possibly update the incumbent
	if(local_lb > status->solution->value) {
        update_best_solution_from_de(status, residual, &desol, local_lb);
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
float get_boole_bound(
        TBKPBBAlgStatus *const status,
        const TBKPBBResidualInstance *const residual,
        size_t current_node,
        size_t* items,
        size_t n_unfixed_items
) {
    TBKPBooleSol boolesol = tbkp_boolesol_lin_gurobi_get(
            status->instance , n_unfixed_items, items, residual->res_capacity);

    float local_lb = (float) (boolesol.lb_sum_profits + residual->sum_profits) *
                     (boolesol.lb_product_probabilities * residual->prod_probabilities);

    if(BB_VERBOSITY_CURRENT > BB_VERBOSITY_INFO) {
        printf("\tLB from Boole relaxation solution: %.3f\n", boolesol.lb);
        printf("\tSum of profits: %" PRIuFAST32 " (fixed %" PRIuFAST32 ")\n",
                boolesol.lb_sum_profits, residual->sum_profits);
        printf("\tProduct of probabilities %f (fixed %f)\n",
                boolesol.lb_product_probabilities, residual->prod_probabilities);
        printf("\tLocal LB: %f (vs %f)\n", local_lb, status->solution->value);
    }

    if(current_node == 1u) {
        status->stats->boole_lb_at_root = local_lb;
        status->stats->time_to_compute_boole_at_root = boolesol.time_to_compute;
    }

    if(local_lb > status->stats->lb) {
        status->stats->lb = local_lb;
    }

    // Possibly update the incumbent
    if(local_lb > status->solution->value) {
        update_best_solution_from_boole(status, residual, &boolesol, local_lb);
    }

    tbkp_boolesol_free_inside(&boolesol);

    return local_lb;
}

/**
 * Solves the deterministic Knapsack Problem when all TB objects are fixed, at a leaf node.
 * If the new solution is better than `status->stats->lb`, it updates it. 
 *
 * @param status                Current state of the algorithm.
 * @param residual              Current residual instance.
 */
static void solve_det_kp(TBKPBBAlgStatus* status, const TBKPBBResidualInstance *const residual) {
		size_t n_det_items = 0u;
		for(size_t i = 0; i < status->instance->n_items; ++i) {
		    if(status->instance->probabilities[i] > 1.0 - EPS) {
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
				if((status->instance->probabilities[i] > 1.0 - EPS)) {
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

			// Compute the new solution value and possibly update the incumbent
			uint_fast32_t sumprof = residual->sum_profits + (uint_fast32_t)myz;
			float zz = (float)sumprof * residual->prod_probabilities;

			if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
			    printf("\tCOMBO solution: %ld, sumprof: %" PRIuFAST32 " => New z: %f (best z: %f)\n",
			            myz, sumprof, zz, status->solution->value);
			}

			if(zz > status->stats->lb) {
                status->stats->lb = zz;
			}

			// Possibly update the incumbent
			if(zz > status->solution->value) {
				status->solution->value = zz;

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
		}
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
        float parent_ub,
        size_t current_node,
        int jbra,
        size_t max_nodes
) {

	// return if the node limit has been reached
    	if ( (max_nodes>0) && ((*(status->n_nodes)) >= max_nodes) ) return;

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

	tbkp_bb_solve_node(status, parent_ub, left_residual, true /* early combo on 1-branch */, max_nodes);

	if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	    printf("[NODE %zu] Returned from left node to father node\n", current_node);
	}

	// return if the node limit has been reached
    	if ( (*(status->n_nodes)) >= max_nodes ) return;

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

	tbkp_bb_solve_node(status, parent_ub, residual, false /* no early combo on 0-branch */, max_nodes);

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
        float parent_ub,
        TBKPBBResidualInstance residual,
        _Bool early_combo,
        size_t max_nodes)
{
    if(timeout(status, parent_ub)) {
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

	float local_ub = parent_ub;
	float local_lb = 0.0f;

	if(status->params->use_de_bounds) {
		DEBounds de_bounds = get_de_bounds(status, &residual, current_node, items, n_unfixed_items);

        if(local_ub == INITIAL_UB_PLACEHOLDER || de_bounds.local_ub < local_ub - EPS) {
            // Local UB improved.
            local_ub = de_bounds.local_ub;
        }

        if(de_bounds.local_lb > local_lb) {
            // Local LB improved.
            local_lb = de_bounds.local_lb;
        }

		if(de_bounds.should_prune) {
			if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
				printf("[NODE %zu] Pruning node thanks to DE bounds.\n", current_node);
			}

			return;
		}
	}

	if(status->params->use_boole_bound && current_node % status->params->boole_bound_frequency == 1u) {
		float boole_lb = get_boole_bound(status, &residual, current_node, items, n_unfixed_items);

		if(boole_lb > local_lb) {
		    // Local LB improved.
		    local_lb = boole_lb;
		}
	}

	if(local_lb > local_ub - EPS && local_lb != INITIAL_LB_PLACEHOLDER && local_ub != INITIAL_UB_PLACEHOLDER) {
	    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
	        printf("[NODE %zu] LB and UB coincide (%.3f). Closing node.\n", current_node, local_lb);
	    }

	    free(items); items = NULL;
	    return;
	}

	if(status->params->use_early_combo && early_combo) {
		// It's smarter to solve the deterministic 01KP early, with the TB items fixed until now. When
		// branching on zero, we don't need to recompute this, because we inherit from the parent not.
		// We only need to call combo again on the 1 branches.

		if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
		    printf("[NODE %zu] Calling COMBO.\n", current_node);
		}

		solve_det_kp(status, &residual);
	}

    // Find the branching item
    int jbra = tbkp_bb_branch_item(status, &residual);

    if(jbra == NO_TIMEBOMB_ITEM_TO_BRANCH) {
        // All uncertain items have been fixed: solve a knapsack instance induced by the unfixed deterministic items.
        // Note: we are in a leaf node!

        if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
            printf("[NODE %zu] No branching item. Leaf node.\n", current_node);
        }

        if(!status->params->use_early_combo) {
            // If not using early combo, we call combo only in the leaf nodes,
            // i.e., here!
            solve_det_kp(status, &residual);
        }

        free(items); items = NULL;
        return;
    }

	branch(status, residual, local_ub, current_node, jbra, max_nodes);

	free(items); items = NULL;
}
/************************************************************
 * END OF LOCAL HELPER FUNCTIONS                            *
 ************************************************************/

TBKPBBSolution* tbkp_branch_and_bound(const TBKPInstance *const instance, TBKPBBStats* stats, TBKPBBParams* params) {
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

    TBKPBBAlgStatus status = {
            .instance = instance,
            .params = params,
            .solution = solution,
            .x = x,
            .stats = stats,
            .n_nodes = &n_nodes
    };

    TBKPBBResidualInstance residual = {
            .prod_probabilities = 1.0f,
            .sum_profits = 0u,
            .res_capacity = instance->capacity
    };

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("Starting solution at the root node.\n");
    }

    tbkp_bb_solve_node(&status, INITIAL_UB_PLACEHOLDER, residual, true /* early combo */, params->num_nodes);

    free(x); x = NULL;

    stats->end_time = clock();
    stats->elapsed_time = (float)(stats->end_time - stats->start_time) / CLOCKS_PER_SEC;
    stats->n_nodes = n_nodes;

    if(BB_VERBOSITY_CURRENT >= BB_VERBOSITY_INFO) {
        printf("Tree exploration completed after %zu nodes.\n", n_nodes);
    }

    if(stats->ub == INITIAL_UB_PLACEHOLDER) {
        // UB was never updated or no open node left

        if(stats->elapsed_time < params->timeout) {
            // Instance solved to optimality
            stats->ub = stats->lb;
        } else {
            // Open instance
            stats->ub = FLT_MAX;
        }
    }

    stats->gap = (stats->ub - stats->lb) / stats->ub;

    return solution;
}
