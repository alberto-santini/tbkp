//
// Created by alberto on 19/02/2020.
//

#include "tbkp_ub.h"
#include "utils/pdqsort_c.h"
#include <stdlib.h>
#include <stdbool.h>
#include <combo.h>
#include <stdio.h>

TBKPDeterministicEqUB tbkp_deub_get(
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

    for(size_t i = 0; i < n_items; ++i) {
        if(cmb_items[i].x) {
            ++n_ub_items;
        }
    }

    size_t* ub_items = malloc(n_ub_items * sizeof(*ub_items));

    for(size_t i = 0; i < n_items; ++i) {
        if((int) (cmb_items[i].x + 0.5)) {
            ub_items[i] = items[cmb_items[i].posizione];
        }
    }

    // TODO: Check if we need the item ids to be sorted.
    // If not, remove the following line.
    pdqsort_c(ub_items, n_items);

    return (TBKPDeterministicEqUB) {.ub = ub, .n_items = n_ub_items, .items = ub_items};
}

void tbkp_deub_free_inside(TBKPDeterministicEqUB *const deub_ptr) {
    free(deub_ptr->items); deub_ptr->items = NULL;
}

void tbkp_deub_print(const TBKPDeterministicEqUB *const ub) {
    printf("UB value: %f\n", ub->ub);
    printf("Objects packed (%zu):\n\t", ub->n_items);
    for(size_t i = 0; i < ub->n_items; ++i) {
        printf("%zu ",ub->items[i]);
    }
    printf("\n");
}
