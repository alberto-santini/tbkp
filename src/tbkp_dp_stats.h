#ifndef TBKP_DP_STATS_H
#define TBKP_DP_STATS_H

#include "tbkp_params.h"
#include <stdbool.h>

typedef struct {
    _Bool built_timebomb_table;
    _Bool built_deterministic_table;
    float elapsed_time;
} TBKPDPStats;

/** Prints the stats to file as a .csv. */
void tbkp_dp_stats_to_file(const TBKPDPStats* stats, const TBKPParams* params);

#endif