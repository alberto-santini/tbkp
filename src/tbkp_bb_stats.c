//
// Created by alberto on 21/04/2020.
//

#include "tbkp_bb_stats.h"
#include "tbkp_bb_defs.h"
#include "tbkp_params.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

TBKPBBStats tbkp_bb_stats_init(void) {
    return (TBKPBBStats) {
            .start_time = clock(),
            .end_time = 0,
            .elapsed_time = 0.0f,
            .lb = INITIAL_LB_PLACEHOLDER,
            .ub = INITIAL_UB_PLACEHOLDER,
            .gap = DBL_MAX,
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

void tbkp_bb_stats_print(const TBKPBBStats *const stats) {
    printf("UB: %.3f; LB: %.3f; Gap: %.3f%%\n", stats->ub, stats->lb, stats->gap * 100.0);
    printf("Elapsed time: %.3f seconds\n", stats->elapsed_time);
    printf("Explored %zu B&B nodes\n", stats->n_nodes);
}

void tbkp_bb_stats_to_file(
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

    if(params->use_cr_bound || params->use_all_bounds_at_root) {
        fprintf(f, ",cr_ub_root,cr_time_root,tot_cr_time,n_cr_called");
    }

    if(params->use_de_bounds || params->use_all_bounds_at_root) {
        fprintf(f, ",de_ub_root,de_lb_root,de_time_root,tot_de_time,n_de_called");
    }

    if(params->use_boole_bound || params->use_all_bounds_at_root) {
        fprintf(f, ",boole_lb_root,boole_time_root,boole_freq,tot_boole_time,n_boole_called");
    }

    fprintf(f, "\n");

    fprintf(f, "%d,%zu,%.3f,%.3f,%.3f,%.3f,%zu",
            params->use_early_combo, params->max_nodes, stats->ub, stats->lb, stats->gap * 100.0, stats->elapsed_time, stats->n_nodes);

    if(params->use_cr_bound || params->use_all_bounds_at_root) {
        fprintf(f, ",%.3f,%.3f", stats->cr_ub_at_root, stats->time_to_compute_cr_at_root);
        fprintf(f, ",%.3f,%zu", stats->tot_time_cr, stats->n_cr_called);
    }

    if(params->use_de_bounds || params->use_all_bounds_at_root) {
        fprintf(f, ",%.3f,%.3f,%.3f", stats->de_ub_at_root, stats->de_lb_at_root, stats->time_to_compute_de_at_root);
        fprintf(f, ",%.3f,%zu", stats->tot_time_de, stats->n_de_called);
    }

    if(params->use_boole_bound || params->use_all_bounds_at_root) {
        fprintf(f, ",%.3f,%.3f,%zu",
                stats->boole_lb_at_root, stats->time_to_compute_boole_at_root, params->boole_bound_frequency);
        fprintf(f, ",%.3f,%zu", stats->tot_time_boole, stats->n_boole_called);
    }

    fprintf(f, "\n");

    fclose(f);
}