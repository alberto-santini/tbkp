#ifndef TBKP_DP_STATS_H
#define TBKP_DP_STATS_H

#include "tbkp_params.h"

typedef struct {
    float elapsed_time;
} TBKPDPStats;

/** Prints the stats to file as a .csv. */
void tbkp_dp_stats_to_file(const TBKPDPStats* stats, const TBKPParams* params);

#endif