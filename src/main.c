#include "tbkp_instance.h"
#include "tbkp_params.h"
#include "tbkp_bb.h"
#include "tbkp_dp.h"
#include <argparse.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gurobi_c.h>

GRBenv* grb_env = NULL;

TBKPParams parse_arguments(int argc, const char** argv) {
    TBKPParams p = {
            .solver = "bb",
            .instance_file = NULL,
            .output_file = NULL,
            .timeout = 3600.0f,
            .boole_bound_frequency = 1u,
            .use_cr_bound = false,
            .use_de_bounds = false,
            .use_boole_bound = false,
            .use_early_combo = false,
            .max_nodes = 0
    };

    // I think we need this because argparse's OPT_BOOLEAN has a bug.
    int early_combo = 0, cont_relax = 0, de_bounds = 0, boole_bound = 0;

    static const char *const usage[] = {
            "tbkp [options]",
            NULL
    };

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_STRING( 'i', "instance", &p.instance_file, "path of the instance"),
        OPT_STRING( 'o', "output", &p.output_file, "path to the csv output file"),
        OPT_STRING( 's', "solver", &p.solver, "solver: bb = branch and bound, dp = dynamic programming"),
        OPT_FLOAT(  't', "timeout", &p.timeout, "timeout in seconds"),
        OPT_INTEGER('c', "earlycombo", &early_combo, "1 if we call COMBO before reaching leaf nodes"),
        OPT_INTEGER('r', "contrelax", &cont_relax, "1 if we use bounds from the continuous relaxation"),
        OPT_INTEGER('d', "debounds", &de_bounds, "1 if we use the DE bounds"),
        OPT_INTEGER('b', "boolebound", &boole_bound, "1 if we use the Boole bound"),
        OPT_INTEGER('f', "boolefreq", &p.boole_bound_frequency, "freequency at which to use the Boole bound"),
        OPT_INTEGER('n', "maxnodes", &p.max_nodes, "Number of branch-and-bound nodes to be explored (1 for root node only, 0 for no limit)"),
        OPT_END()
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argc = argparse_parse(&argparse, argc, argv);

    p.use_early_combo = (early_combo == 1);
    p.use_cr_bound = (cont_relax == 1);
    p.use_de_bounds = (de_bounds == 1);
    p.use_boole_bound = (boole_bound == 1);

    if(!p.instance_file) {
        printf("You have to specify an instance file!\n");
        exit(EXIT_FAILURE);
    }

    return p;
}

int main(int argc, const char** argv) {
    TBKPParams p = parse_arguments(argc, argv);
    TBKPInstance* instance = tbkp_instance_read(p.instance_file);

    if(strcmp(p.solver, "bb") == 0) {
        TBKPBBStats stats = tbkp_stats_init();
        int error = GRBloadenv(&grb_env, NULL);

        if(!p.use_de_bounds) {
            printf("Error: must use DE bounds.\n");
            exit(EXIT_FAILURE);
        }

        if(error) {
            printf("Gurobi loadenv error: %d\n", error);
            exit(EXIT_FAILURE);
        }

        error = GRBsetintparam(grb_env, GRB_INT_PAR_OUTPUTFLAG, 0);

        if(error) {
            printf("Gurobi setintparam output flag error: %d\n", error);
            exit(EXIT_FAILURE);
        }
        
        TBKPBBSolution* bbsol = tbkp_branch_and_bound(instance, &stats, &p);
        tbkp_stats_to_file(&stats, &p);

        // Clean up
        tbkp_sol_free(&bbsol);
        GRBfreeenv(grb_env);
    } else if(strcmp(p.solver, "dp") == 0) {
        const float solution = tbkp_dp_solve(instance);

        printf("Dynamic programming solution: %.3f\n", solution);
    } else {
        printf("Solver not recognised: %s\n", p.solver);
        
        tbkp_instance_free(&instance); instance = NULL;
        exit(EXIT_FAILURE);
    }

    tbkp_instance_free(&instance); instance = NULL;

    return 0;
}
