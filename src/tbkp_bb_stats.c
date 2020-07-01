//
// Created by alberto on 21/04/2020.
//

#include "tbkp_bb_stats.h"
#include "tbkp_bb_defs.h"
#include "tbkp_params.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

TBKPBBStats tbkp_stats_init(void) {
    return (TBKPBBStats) {
            .start_time = clock(),
            .end_time = 0,
            .elapsed_time = 0.0f,
            .lb = INITIAL_LB_PLACEHOLDER,
            .ub = INITIAL_UB_PLACEHOLDER,
            .gap = FLT_MAX,
            .boole_lb_at_root = INITIAL_LB_PLACEHOLDER,
            .de_lb_at_root = INITIAL_LB_PLACEHOLDER,
            .de_ub_at_root = INITIAL_UB_PLACEHOLDER,
            .cr_ub_at_root = INITIAL_UB_PLACEHOLDER,
            .time_to_compute_cr_at_root = 0.0f,
            .time_to_compute_boole_at_root = 0.0f,
            .time_to_compute_de_at_root = 0.0f,
            .n_nodes = 0u,
    };
}

void tbkp_stats_print(const TBKPBBStats *const stats) {
    printf("UB: %.3f; LB: %.3f; Gap: %.3f%%\n", stats->ub, stats->lb, stats->gap * 100.0f);
    printf("Elapsed time: %.3f seconds\n", stats->elapsed_time);
    printf("Explored %zu B&B nodes\n", stats->n_nodes);
}

void tbkp_stats_to_file(
        const TBKPBBStats *const stats,
        const TBKPParams *const params)
{
    if(!params->output_file) {
        return;
    }

    FILE* f = fopen(params->output_file, "w");

    if(!f) {
        printf("Error opening csv output file %s\n", params->output_file);
        exit(EXIT_FAILURE);
    }

    fprintf(f, "early_combo,max_nodes,ub,lb,gap,time_s,n_nodes");

    if(params->use_cr_bound) {
        fprintf(f, ",cr_ub_root,cr_time_root");
    }

    if(params->use_de_bounds) {
        fprintf(f, ",de_ub_root,de_lb_root,de_time_root");
    }

    if(params->use_boole_bound) {
        fprintf(f, ",boole_lb_root,boole_time_root,boole_freq");
    }

    fprintf(f, "\n");

    fprintf(f, "%d,%zu,%.3f,%.3f,%.3f,%.3f,%zu",
            params->use_early_combo, params->max_nodes, stats->ub, stats->lb, stats->gap * 100.0f, stats->elapsed_time, stats->n_nodes);

    if(params->use_cr_bound) {
        fprintf(f, ",%.3f,%.3f", stats->cr_ub_at_root, stats->time_to_compute_cr_at_root);
    }

    if(params->use_de_bounds) {
        fprintf(f, ",%.3f,%.3f,%.3f", stats->de_ub_at_root, stats->de_lb_at_root, stats->time_to_compute_de_at_root);
    }

    if(params->use_boole_bound) {
        fprintf(f, ",%.3f,%.3f,%zu",
                stats->boole_lb_at_root, stats->time_to_compute_boole_at_root, params->boole_bound_frequency);
    }

    fprintf(f, "\n");

    fclose(f);
}