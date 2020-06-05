//
// Created by alberto on 21/04/2020.
//

#ifndef TBKP_TBKP_BB_STATS_H
#define TBKP_TBKP_BB_STATS_H

#include "tbkp_bb_params.h"
#include <stddef.h>
#include <stdbool.h>
#include <time.h>

/** Collects statistics on the B&B solution process. */
typedef struct {
    clock_t start_time;
    clock_t end_time;
    float elapsed_time;
    float ub; /* Starts at value -1.0f */
    float lb;
    float gap;
    float boole_lb_at_root;
    float de_lb_at_root;
    float de_ub_at_root;
    float time_to_compute_boole_at_root;
    float time_to_compute_de_at_root;
    size_t n_nodes;
} TBKPBBStats;

/** Initialises an empty statistics object. */
TBKPBBStats tbkp_stats_init(void);

/** Prints a summary of the statistics to stdout. */
void tbkp_stats_print(const TBKPBBStats* stats);

/** Prints the stats to file as a .csv. */
void tbkp_stats_to_file(const TBKPBBStats* stats, const TBKPBBParams* params);

#endif //TBKP_TBKP_BB_STATS_H