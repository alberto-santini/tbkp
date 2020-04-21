//
// Created by alberto on 21/04/2020.
//

#ifndef TBKP_TBKP_BB_STATS_H
#define TBKP_TBKP_BB_STATS_H

#include <stddef.h>
#include <stdbool.h>
#include <time.h>

/** Collects statistics on the B&B solution process. */
typedef struct {
    clock_t start_time;
    clock_t end_time;
    float timeout;
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
    size_t boole_bound_frequency;
    _Bool use_de_bounds;
    _Bool use_boole_bound;
    _Bool use_early_combo;
} TBKPBBStats;

/** Initialises an empty statistics object. */
TBKPBBStats tbkp_stats_init(
        float timeout,
        _Bool use_early_combo,
        _Bool use_de_bounds,
        _Bool use_bool_bound,
        size_t boole_bound_frequency);

/** Prints a summary of the statistics to stdout. */
void tbkp_stats_print(const TBKPBBStats* stats);

/** Prints the stats to file as a .csv. */
void tbkp_stats_to_file(const TBKPBBStats* stats, const char* csv_filename);

#endif //TBKP_TBKP_BB_STATS_H
