//
// last update: Michele 03/03/2020.
//

#include "tbkp_BaB.h"
#include "utils/pdqsort_c.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <combo.h>
#include "tbkp_de_sol.h"

#define PRILEV 100
#define EPS 1e-6

int branch_and_bound(const TBKPInstance *const instance)
{
	TBKPsolution* solution = malloc(sizeof(*solution));
        solution->x = malloc(instance->n_items * sizeof(int));
        solution->prod_Probab = 1.0;
        solution->sum_Profits = 0;
        solution->value = 0.0;

	// x[j] >= 0 if item j has been fixed by branching
	//      = -1 if item j is still unfixed
	int* x = malloc(instance->n_items * sizeof(*x));
	for ( size_t t = 0; t < instance->n_items; ++t) x[t] = -1;
	float prod_Probab = 1.0;
	uint_fast32_t sum_Profits = 0;
	uint_fast32_t resCapa = instance->capacity;

	int nnodes = 0;

	if ( PRILEV >= 1000 )
	{
		printf("root node; chiamo solve_node con capacity %d e n_items %d\n", (int) instance->capacity, (int) instance->n_items);
		for ( size_t t = 0; t < instance->n_items; t++ ) printf(" %3d: prof %4d, weight %4d, probab %7.3lf\n", (int) t, instance->profits[t], instance->weights[t], instance->probabilities[t]);
		printf("\n");
	}

	solve_node(instance, &nnodes, x, prod_Probab, sum_Profits, resCapa, solution);

	printf("\nbranch-and-bound: solution value %f\n", solution->value);
	for ( size_t t = 0; t < instance->n_items; t++ ) printf(" %3d: prof %4d, weight %4d, probab %7.3lf --> x %1d\n", (int) t, instance->profits[t], instance->weights[t], instance->probabilities[t], solution->x[t]);

	return 1;
}



// returns the index of the item selected for branching (-1 if NONE)
int branch_item(const TBKPInstance *const instance, int *x, uint_fast32_t capacity)
{
	// find the first uncertain item that is not fixed and fits in the residual capacity
	for ( size_t i = 0; i < instance->n_items; ++i)
	{
		if ( instance->probabilities[i] > 1.0 - EPS ) continue;
		if ( x[i] >= 0 ) continue;
		if ( instance->weights[i] > capacity ) continue;
		return i;
	}
	return -1;
}


void solve_node(const TBKPInstance *const instance, int *nnodes, int* x, float prod_Probab, uint_fast32_t sum_Profits, uint_fast32_t resCapa, TBKPsolution *solution)
{
	(*nnodes)++;
	int currnode = *nnodes;

	// count the number of items that are still unfixed
	size_t n_items = 0;
	for ( size_t j = 0; j < instance->n_items; j++ ) if ( x[j] < 0 ) n_items++;
	if ( n_items <= 0 ) return;

	// define the residual instance
	size_t* items = malloc(n_items * sizeof(*items));
	int cnt = 0;
	for ( size_t j = 0; j < instance->n_items; j++ ) 
	{
		if ( x[j] < 0 ) 
		{
			items[cnt] = j;
			cnt++;
		}
	}	
	
	if ( PRILEV >= 1000 ) 
	{
		printf("\n\nsono nel nodo %d con bestZ %f e resCapa %d\n", currnode, solution->value, (int) resCapa);
		printf("sumProfits %d, productProfits %f\n", (int) sum_Profits, prod_Probab);
		printf("number of residual items %d: ", (int) n_items);
		for ( size_t j = 0; j < n_items; j++ ) printf("%2d ", (int) items[j]);
		printf("\n");
	}

	// solve the deterministic relaxation
	TBKPDeterministicEqSol desol = tbkp_desol_get(instance, n_items, items, resCapa);
	if ( PRILEV >= 1000 ) tbkp_desol_print(&desol);

	// compute the local upper bound and possibly kill the node
	float localub = (sum_Profits + desol.ub) * prod_Probab;
	if ( PRILEV >= 1000 ) printf("desol.ub %f, sumProfits %d, productProfits %f -> localub %f (vs %f)\n", desol.ub, (int) sum_Profits, prod_Probab, localub, solution->value);

	if ( localub <= solution->value ) 
	{
		if ( PRILEV >= 1000 ) printf("node killed\n");
		free(items);
		return;
	}



	// possibly update the incumbent
	if ( desol.lb > solution->value ) {
		solution->value = desol.lb;
		for ( size_t i = 0; i < instance->n_items; ++i) 
		for ( size_t i = 0; i < instance->n_items; ++i) 
		{
			if ( x[i] >= 0 ) solution->x[i] = x[i];
		}
		for ( size_t i = 0; i < desol.n_items; ++i)
		{
			size_t j = desol.items[i];
			solution->x[j] = 1;
		}
		printf("solution update: new value %f\n", solution->value);
	}

	// find the branching item
	int jbra = branch_item(instance, x, resCapa);
	if ( jbra < 0 ) 
	{
		if ( PRILEV >= 1000 ) printf("no branching item\n");

		// all uncertain items have been fixed: solve a knapsack instance induced by the unfixed deterministic items
		int n_items = 0;
		for ( size_t i = 0; i < instance->n_items; ++i) if ( (instance->probabilities[i] > 1.0 - EPS) && (x[i] < 0) ) n_items++;
		if ( n_items > 0 ) 
		{
			cmb_item* cmb_items = malloc(n_items * sizeof(*cmb_items));
			cmb_stype sumW = 0;
			cmb_stype sumP = 0;
			int cnt = 0;
			for ( size_t i = 0; i < instance->n_items; ++i ) 
			{
				if ( (instance->probabilities[i] > 1.0 - EPS) && (x[i] < 0) ) 
				{
					cmb_items[cnt] = (cmb_item) 
					{
						.p = (cmb_itype) instance->profits[i],
						.w = (cmb_itype) instance->weights[i],
						.x = 0,
						.pos = (int) i
					};
					sumW += cmb_items[cnt].w;
					sumP += cmb_items[cnt].p;
					cnt++;
				}
			}
			int myz = 0;
			if ( sumW > resCapa ) 
			{
				const long cmb_z = combo(&cmb_items[0], &cmb_items[n_items - 1], (cmb_stype)resCapa, 0, INT32_MAX, true, false);
				myz = cmb_z;
			}
			else
			{	
				myz = sumP;
				for (int cnt = 0; cnt < n_items; ++cnt ) cmb_items[cnt].x = 1;
			}

			// compute the new solution value and possibly update the incumbent
			const long sumprof = sum_Profits + myz;
			double zz = sumprof * prod_Probab;	

			if ( PRILEV >= 1000 ) printf("myz %d, sumprof %d -> zz %f (bestZ %f)\n", (int) myz, (int) sumprof, zz, solution->value);

			// possibly update the incumbent
			if ( zz > solution->value ) 
			{
				solution->value = zz;
				for ( size_t i = 0; i < instance->n_items; i++ ) 
				{
					if ( x[i] == 1 ) solution->x[i] = 1;
					else solution->x[i] = 0;
				}
				for (int cnt = 0; cnt < n_items; ++cnt ) 
				{
					if ( cmb_items[cnt].x > 0.5 ) 
					{
						int i = cmb_items[cnt].pos;
						solution->x[i] = 1;
					}
				}
			}
		}
		free(items);
		return;
	}

	// Left node: fix the item in the solution
	uint_fast32_t new_resCapa = resCapa - instance->weights[jbra];
	uint_fast32_t newsum_Profits = sum_Profits + instance->profits[jbra];
	float newProd_Probab = prod_Probab * instance->probabilities[jbra];
	x[jbra] = 1;
	if ( PRILEV >= 1000 ) 
	{
		printf("node %d; item jbra=%d (prof %d, peso %d, prob %f) fissato\n", currnode, (int) jbra, (int) instance->profits[jbra], (int) instance->weights[jbra], instance->probabilities[jbra]);
		printf("chiamo solve_node con new_resCapa %d, ", (int) new_resCapa);
		printf("newsum_Profits %d, newProd_Probab %f\n", (int) newsum_Profits, newProd_Probab);
	}

	solve_node(instance, nnodes, x, newProd_Probab, newsum_Profits, new_resCapa, solution);

	if ( PRILEV >= 1000 ) printf("\n\ntornato dal left node: nodo %d (nnodes %d)\n", currnode, *nnodes);

	// Right node: remove the item
	new_resCapa = resCapa;
	newsum_Profits = sum_Profits;
	newProd_Probab = prod_Probab;
	x[jbra] = 0;
	if ( PRILEV >= 1000 ) 
	{
		printf("node %d; item jbra=%d (prof %d, peso %d, prob %f) rimosso\n", currnode, (int) jbra, (int) instance->profits[jbra], (int) instance->weights[jbra], instance->probabilities[jbra]);
		printf("chiamo solve_node con new_resCapa %d, ", (int) new_resCapa);
		printf("newsum_Profits %d, newProd_Probab %f\n", (int) newsum_Profits, newProd_Probab);
	}

	solve_node(instance, nnodes, x, newProd_Probab, newsum_Profits, new_resCapa, solution);

	if ( PRILEV >= 1000 ) printf("\n\ntornato dal right node: nodo %d (nnodes %d)\n", currnode, *nnodes);

	x[jbra] = -1;

	free(items);
	return;
}
