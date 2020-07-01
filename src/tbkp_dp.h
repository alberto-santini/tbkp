#ifndef TBKP_DP_H
#define TBKP_DP_H

#include <stdlib.h>
#include "tbkp_instance.h"

#define TB_TABLE_UNKNOWN -1.0f
#define DET_TABLE_UNKNOWN UINT_FAST32_MAX

typedef struct {
    float* tb_table;
    uint_fast32_t* det_table;
} TBKPDPTables;

TBKPDPTables tbkp_dp_tables_init(const TBKPInstance* inst);
void tbkp_dp_tables_free(TBKPDPTables* t);

#endif