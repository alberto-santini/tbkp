//
// Created by alberto on 23/04/2020.
//

#ifndef TBKP_TBKP_PARAMS_H
#define TBKP_TBKP_PARAMS_H

#include <stddef.h>

typedef struct {
    char* solver;
    char* instance_file;
    char* output_file;
    float timeout;

    // BB options:
    size_t boole_bound_frequency;
    _Bool use_cr_bound;
    _Bool use_de_bounds;
    _Bool use_boole_bound;
    _Bool use_early_combo;
    _Bool use_early_pruning;
    _Bool use_all_bounds_at_root;
    size_t max_nodes;
    float boole_lin_solver_timeout_s;
} TBKPParams;

#endif //TBKP_TBKP_PARAMS_H
