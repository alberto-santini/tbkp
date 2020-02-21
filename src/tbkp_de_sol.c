//
// Created by alberto on 19/02/2020.
//

#include "tbkp_de_sol.h"
#include "utils/pdqsort_c.h"
#include <stdlib.h>
#include <stdbool.h>
#include <combo.h>
#include <stdio.h>

TBKPDeterministicEqSol tbkp_desol_get(
        const TBKPInstance *const instance, size_t n_items,
        const size_t *const items, uint_fast32_t capacity)
{
    // Multiplier we need because COMBO only takes integer profits!
    const float cmb_multiplier = 10000.0f;

    // AS: Michele, è corretto (e consigliato). Guarda qua:
    // https://stackoverflow.com/questions/17258647/why-is-it-safer-to-use-sizeofpointer-in-malloc
    item* cmb_items = malloc(n_items * sizeof(*cmb_items));

    for(size_t i = 0; i < n_items; ++i) {
        const float real_w = (float)(instance->profits[items[i]]) * instance->probabilities[items[i]];
        cmb_items[i] = (item) {
            .p = (itype) (real_w * cmb_multiplier),
            .w = (itype) instance->weights[items[i]],
            .x = 0,
            .posizione = (int) i
        };
    }

    const long cmb_ub = combo(&cmb_items[0], &cmb_items[n_items - 1], (stype)capacity, 0, INT32_MAX, true, false);
    const float ub = (float)cmb_ub / cmb_multiplier;

    size_t n_ub_items = 0;
    float lb = 0.0f;

    // We use this loop to count how many items are packed by COMBO, but also to compute
    // the first part of the 01-KP objective function (the sum of the profits).
    for(size_t i = 0; i < n_items; ++i) {
        if(cmb_items[i].x) {
            ++n_ub_items;
            lb += (float) instance->profits[items[cmb_items[i].posizione]];
        }
    }

    size_t* ub_items = malloc(n_ub_items * sizeof(*ub_items));

    // We use this loop to create the list of the items packed by COMBO, but also to compute
    // the second part of the 01-KP objective function (the product of the probabilities).
    for(size_t i = 0; i < n_items; ++i) {
        if(cmb_items[i].x) {
            ub_items[i] = items[cmb_items[i].posizione];
            lb *= instance->probabilities[items[cmb_items[i].posizione]];
        }
    }

    // TODO: Check if we need the item ids to be sorted.
    // If not, remove the following line.
    pdqsort_c(ub_items, n_items);

    return (TBKPDeterministicEqSol) {.ub = ub, .lb = lb, .n_items = n_ub_items, .items = ub_items};
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
