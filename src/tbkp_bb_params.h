//
// Created by alberto on 23/04/2020.
//

#ifndef TBKP_TBKP_BB_PARAMS_H
#define TBKP_TBKP_BB_PARAMS_H

#include <stddef.h>

typedef struct {
    char* instance_file;
    char* output_file;
    float timeout;
    size_t boole_bound_frequency;
    _Bool use_de_bounds;
    _Bool use_boole_bound;
    _Bool use_early_combo;
} TBKPBBParams;

#endif //TBKP_TBKP_BB_PARAMS_H
