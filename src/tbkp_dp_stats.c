#include "tbkp_dp_stats.h"
#include "tbkp_params.h"
#include <stdio.h>
#include <stdlib.h>

void tbkp_dp_stats_to_file(const TBKPDPStats* stats, const TBKPParams* params) {
    if(!params->output_file) {
        return;
    }

    FILE* f = fopen(params->output_file, "w");

    if(!f) {
        printf("Error opening csv output file %s\n", params->output_file);
        exit(EXIT_FAILURE);
    }

    fprintf(f, "time_s,built_tb,built_det\n%.2f,%d,%d\n",
            stats->elapsed_time, stats->built_timebomb_table, stats->built_deterministic_table);

    fclose(f);
}