//
// Created by alberto on 18/02/2020.
//

#include "tbkp_instance.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <combo.h>
#include <assert.h>

#define EPS 1e-6

typedef struct {
    double pi;
    uint_fast32_t w;
    uint_fast32_t p;
} Item;

static int compare_items_by_probs(const void* x, const void* y) {
    const Item* it_x = x;
    const Item* it_y = y;
    return it_x->pi < it_y->pi ? -1 : it_x->pi > it_y->pi;
}

static uint_fast32_t tb_items_ub(const TBKPInstance *const inst) {
    cmb_item* cmb_items = malloc(inst->n_tb_items * sizeof(*cmb_items));
    cmb_stype tot_p = 0;
    
    for(size_t i = 0u; i < inst->n_tb_items; ++i) {
        assert(inst->probabilities[i] < 1.0 - EPS);
        cmb_items[i] = (cmb_item){
            .p = (cmb_itype) inst->profits[i],
            .w = (cmb_itype) inst->weights[i],
            .x = 0,
            .pos = i
        };
        tot_p += (cmb_stype)inst->profits[i];
    }

    if(inst->tb_tot_weight <= inst->capacity) {
        free(cmb_items);
        return (uint_fast32_t)tot_p;
    }

    uint_fast32_t sol = (uint_fast32_t) combo(
        &cmb_items[0], &cmb_items[inst->n_tb_items - 1],
        (cmb_stype)inst->capacity, 0, INT32_MAX, true, false);

    free(cmb_items);

    return sol;
}

/** Arrange items by increasing probabilities. As a side-effect, all time-bomb
 *  items are at the beginning of the lists. It also sets last_tb_item_index.
 */
static void sort_instance_by_increasing_probabilities(TBKPInstance* inst, Item* items) {
    qsort(items, inst->n_items, sizeof(*items), compare_items_by_probs);

    for(size_t i = 0u; i < inst->n_items; ++i) {
        inst->weights[i] = items[i].w;
        inst->profits[i] = items[i].p;
        inst->probabilities[i] = items[i].pi;

        if(items[i].pi < 1.0) {
            inst->last_tb_item_index = i;
        }
    }
}

TBKPInstance* tbkp_instance_read(const char *const filename) {
    TBKPInstance* instance = malloc(sizeof(*instance));

    if(instance == NULL) {
        printf("Cannot allocate memory for the instance\n");
        exit(EXIT_FAILURE);
    }

    FILE* fd;
    fd = fopen(filename, "r");

    if(fd == NULL) {
        printf("Cannot read %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int n_read;

    n_read = fscanf(fd, "%" SCNuFAST32 " %" SCNuFAST32, &(instance->n_items), &(instance->capacity));

    if(n_read != 2) {
        printf("Error reading the first line with number of items and capacity\n");
        exit(EXIT_FAILURE);
    }

    if(instance->n_items == 0) {
        printf("Error: you have an instance with no items\n");
        exit(EXIT_FAILURE);
    }

    if(instance->n_items > instance->capacity) {
        printf("Error: you have an instance with more items than capacity units\n");
        exit(EXIT_FAILURE);
    }

    Item* items = malloc(instance->n_items * sizeof(*items));
    instance->weights = malloc(instance->n_items * sizeof(*(instance->weights)));
    instance->profits = malloc(instance->n_items * sizeof(*(instance->profits)));
    instance->probabilities = malloc(instance->n_items * sizeof(*(instance->probabilities)));
    instance->n_tb_items = 0u;
    instance->n_det_items = 0u;
    instance->tb_tot_weight = 0u;
    instance->tb_ub_packed_profit = 0u;
    instance->sum_profit_times_probability = 0.0;

    if((instance->weights == NULL) || (instance->profits == NULL) || (instance->probabilities == NULL)) {
        printf("Cannot allocate memory for data\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < instance->n_items; ++i) {
        n_read = fscanf(fd, "%" SCNuFAST32 " %" SCNuFAST32 " %lf",
                &(items[i].w), &(items[i].p), &(items[i].pi));

        if(n_read != 3) {
            printf("Error reading item %zu\n", i);
            exit(EXIT_FAILURE);
        }

        if(items[i].w > instance->capacity) {
            printf("Wrong weight for item %zu: %" PRIuFAST32 " larger than capacity %" PRIuFAST32 "\n",
                    i, instance->weights[i], instance->capacity);
            exit(EXIT_FAILURE);
        }

        if(items[i].pi < 0.0 || items[i].pi > 1.0) {
            printf("Wrong probability for item %zu: %f\n", i, instance->probabilities[i]);
            exit(EXIT_FAILURE);
        }

        if(items[i].pi < 1.0 - EPS) {
            ++(instance->n_tb_items);
            instance->tb_tot_weight += items[i].w;
        } else {
            ++(instance->n_det_items);
            instance->det_tot_weight += items[i].w;
        }

        instance->sum_profit_times_probability += (double) items[i].p * items[i].pi;
    }

    if(fclose(fd) != 0) {
        printf("Error closing file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    sort_instance_by_increasing_probabilities(instance, items);
    instance->tb_ub_packed_profit = tb_items_ub(instance);

    free(items);

    return instance;
}

void tbkp_instance_free(TBKPInstance** instance_ptr) {
    free((*instance_ptr)->weights); (*instance_ptr)->weights = NULL;
    free((*instance_ptr)->profits); (*instance_ptr)->profits = NULL;
    free((*instance_ptr)->probabilities); (*instance_ptr)->probabilities = NULL;
    free(*instance_ptr); *instance_ptr = NULL;
}

void tbkp_instance_print(const TBKPInstance *const instance) {
    printf("Num items: %" PRIuFAST32 ", Capacity: %" PRIuFAST32 "\n", instance->n_items, instance->capacity);
    printf("Last TB item: %" PRIuFAST32 "\n", instance->last_tb_item_index);

    for(size_t i = 0; i < instance->n_items; ++i) {
        printf("Object %zu: weight=%" PRIuFAST32 ", profit=%" PRIuFAST32 ", prob=%f\n",
                i, instance->weights[i], instance->profits[i], instance->probabilities[i]);
    }
}