//
// Created by alberto on 19/02/2020.
// Modified by Michele on 3/3/2020
//

#include "tbkp_de_sol.h"
#include "utils/pdqsort_c.h"
#include <stdlib.h>
#include <stdbool.h>
#include <combo.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <assert.h>

TBKPDeterministicEqSol tbkp_desol_get(
        const TBKPInstance *const instance, size_t n_items,
        const size_t *const items, uint_fast32_t capacity)
{
    clock_t start_time = clock();

    // Multiplier we need because COMBO only takes integer profits!
    const double cmb_multiplier = 1000.0;

    cmb_item* cmb_items = malloc(n_items * sizeof(*cmb_items));

    if(!cmb_items) {
        printf("Cannot allocate memory for combo items in tbkp_desol_get!\n");
        exit(EXIT_FAILURE);
    }

    cmb_stype sumw = 0;
    cmb_stype sump = 0;
    unsigned long long total_profit = 0ul;

    for(size_t i = 0; i < n_items; ++i) {
        const double real_w = (double)(instance->profits[items[i]]) * instance->probabilities[items[i]];
        cmb_items[i] = (cmb_item) {
            .p = (cmb_itype) (real_w * cmb_multiplier),
            .w = (cmb_itype) instance->weights[items[i]],
            .x = 0,
            .pos = i
        };

        assert(cmb_items[i].p >= 0);
        total_profit += (unsigned long long) cmb_items[i].p;

        sump += cmb_items[i].p;
        sumw += cmb_items[i].w;
    }

    if(total_profit > LONG_MAX) {
        printf("Total profits (including multipliers) too large! Increase the size of int used by combo.\n");
        exit(EXIT_FAILURE);
    }

    cmb_stype cmb_ub;
    
    // Combo crashes if the total weights are < than the knapsack's capacity!
    if(sumw <= (cmb_stype)capacity) {
        cmb_ub = sump;
        for(size_t i = 0; i < n_items; ++i) {
            cmb_items[i].x = 1;
        }
    } else {
        cmb_ub = combo(&cmb_items[0], &cmb_items[n_items - 1], (cmb_stype)capacity, 0, INT32_MAX, true, false);
    }

    const double ub = (double)cmb_ub / cmb_multiplier;

    size_t n_ub_items = 0;
    uint_fast32_t lb_sum = 0;

    // We use this loop to count how many items are packed by COMBO, but also to compute
    // the first part of the 01-KP objective function (the sum of the profits).
    for(size_t i = 0; i < n_items; ++i) {
        if(cmb_items[i].x) {
            ++n_ub_items;
            lb_sum += instance->profits[items[cmb_items[i].pos]];
        }
    }

    size_t* ub_items = malloc(n_ub_items * sizeof(*ub_items));

    if(!ub_items) {
        printf("Cannot allocate memory for UB items\n");
        exit(EXIT_FAILURE);
    }

    // We use this loop to create the list of the items packed by COMBO, but also to compute
    // the second part of the 01-KP objective function (the product of the probabilities).
    size_t curr_id = 0u;
    double lb_prod = 1.0;
    for(size_t i = 0; i < n_items; ++i) {
        if(cmb_items[i].x) {
            ub_items[curr_id++] = items[cmb_items[i].pos];
            lb_prod *= instance->probabilities[items[cmb_items[i].pos]];
        }
    }

    // TODO: Check if we need the item ids to be sorted.
    // If so, uncomment the following line.
    // pdqsort_c(ub_items, n_ub_items);

    free(cmb_items); cmb_items = NULL;

    clock_t end_time = clock();
    float elapsed_time = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    return (TBKPDeterministicEqSol) {
        .ub = ub,
        .lb = (double) lb_sum * lb_prod,
        .lb_sum_profits = lb_sum,
        .lb_product_probabilities = lb_prod,
        .n_items = n_ub_items,
        .items = ub_items,
        .time_to_compute = elapsed_time
    };
}

void tbkp_desol_free_inside(TBKPDeterministicEqSol *const desol_ptr) {
    free(desol_ptr->items); desol_ptr->items = NULL;
}

void tbkp_desol_print(const TBKPDeterministicEqSol *const sol) {
    printf("UB value: %f\n", sol->ub);
    printf("LB value: %f\n", sol->lb);
    printf("Objects packed (%zu):\n\t", sol->n_items);
    for(size_t i = 0; i < sol->n_items; ++i) {
        printf("%zu ", sol->items[i]);
    }
    printf("\n");
}
