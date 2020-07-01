//
// Created by alberto on 18/02/2020.
//

#ifndef TBKP_TBKP_INSTANCE_H
#define TBKP_TBKP_INSTANCE_H

#include <inttypes.h>

typedef struct {
    uint_fast32_t n_items;
    uint_fast32_t capacity;
    uint_fast32_t* weights;
    uint_fast32_t* profits;
    float* probabilities;

    uint_fast32_t n_tb_items;
    uint_fast32_t n_det_items;
    uint_fast32_t last_tb_item_index;

    uint_fast32_t tb_tot_weight;
    uint_fast32_t tb_ub_packed_profit;
    uint_fast32_t det_tot_weight;
} TBKPInstance;

TBKPInstance* tbkp_instance_read(const char* filename);
void tbkp_instance_free(TBKPInstance** instance_ptr);
void tbkp_instance_print(const TBKPInstance* instance);

#endif //TBKP_TBKP_INSTANCE_H
