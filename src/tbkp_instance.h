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
} TBKPInstance;

TBKPInstance* tbkp_instance_read(const char* filename);
void tbkp_instance_free(TBKPInstance** instance_ptr);
void tbkp_instance_print(const TBKPInstance* instance);

#endif //TBKP_TBKP_INSTANCE_H
