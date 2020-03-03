//
// Created by alberto on 19/02/2020.
//

#ifndef TBKP_TBKP_UB_H
#define TBKP_TBKP_UB_H

#include <stddef.h>
#include "tbkp_instance.h"

typedef struct {
    int* x;
    float prod_Probab;
    uint_fast32_t sum_Profits;
    float value;
} TBKPsolution;

int branch_and_bound(const TBKPInstance *const instance);
int branch_item(const TBKPInstance *const instance, int *x, uint_fast32_t capacity);

uint_fast32_t compute_LB(const TBKPInstance *const instance, size_t n_items, const size_t *const items, uint_fast32_t capacity);
void solve_node(const TBKPInstance *const instance, int *nnodes, int* x, float prod_Probab, uint_fast32_t sum_Profit, uint_fast32_t resCapa, TBKPsolution *solution);

#endif //TBKP_TBKP_UB_H
