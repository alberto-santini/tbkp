//
// Created by alberto on 21/04/2020.
//

#ifndef TBKP_TBKP_BB_STATS_H
#define TBKP_TBKP_BB_STATS_H

#include "tbkp_params.h"
#include "tbkp_bb_solution.h"
#include <stddef.h>
#include <stdbool.h>
#include <time.h>

/** Collects statistics on the B&B solution process. */
typedef struct {
    clock_t start_time;
    clock_t end_time;
    float elapsed_time;

    double ub; /* Starts at value -1.0 */
    double lb;
    double gap;

    double boole_lb_at_root;
    double de_lb_at_root;
    double de_ub_at_root;
    double cr_ub_at_root;

    float time_to_compute_boole_at_root;
    float time_to_compute_de_at_root;
    float time_to_compute_cr_at_root;

    float tot_time_boole;
    float tot_time_de;
    float tot_time_cr;

    size_t n_boole_called;
    size_t n_de_called;
    size_t n_cr_called;

    size_t n_nodes;
} TBKPBBStats;

/** Initialises an empty statistics object. */
TBKPBBStats tbkp_bb_stats_init(void);

/** Prints a summary of the statistics to stdout. */
void tbkp_bb_stats_print(const TBKPBBStats* stats);

/** Prints the stats to file as a .csv. */
void tbkp_bb_stats_to_file(const TBKPBBStats* stats, const TBKPBBSolution* sol, const TBKPInstance* instance, const TBKPParams* params);

#endif //TBKP_TBKP_BB_STATS_H
