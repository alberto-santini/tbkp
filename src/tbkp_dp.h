#ifndef TBKP_DP_H
#define TBKP_DP_H

#include <stdlib.h>
#include "tbkp_dp_stats.h"
#include "tbkp_instance.h"

typedef struct {
    float* rows[2];
    size_t active_row;
    uint_fast32_t last_index;
} TBKPDPTimeBombTable;

typedef struct {
    uint_fast32_t* t;
} TBKPDPDeterministicTable;

TBKPDPTimeBombTable tbkp_dp_timebombtable_init(const TBKPInstance* inst);
void tbkp_dp_timebombtable_free(TBKPDPTimeBombTable* tb_t);

TBKPDPDeterministicTable tbkp_dp_deterministictable_init(const TBKPInstance* inst);
void tbkp_dp_deterministictable_free(TBKPDPDeterministicTable* d_t);

float tbkp_dp_solve(const TBKPInstance *const inst, TBKPDPStats* stats, const TBKPParams *const params);

#endif