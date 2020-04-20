#include "tbkp_instance.h"
#include "tbkp_de_sol.h"
#include "tbkp_boole_sol.h"
#include "tbkp_BaB.h"
#include <argparse.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <gurobi_c.h>

GRBenv* grb_env = NULL;

typedef struct Params {
    char* instance_file;
    char* output_file;
    float timeout_s;
    size_t boole_bound_freq;
    _Bool de_bounds;
    _Bool boole_bound;
} Params;

Params parse_arguments(int argc, const char** argv) {
    Params p = {
            .instance_file = NULL,
            .output_file = NULL,
            .timeout_s = 3600.0f,
            .boole_bound_freq = 1u,
            .de_bounds = false,
            .boole_bound = false
    };

    static const char *const usage[] = {
            "tbkp [options]",
            NULL
    };

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_STRING('i', "instance", &p.instance_file, "path of the instance"),
        OPT_STRING('o', "output", &p.output_file, "path to the csv output file"),
        OPT_FLOAT('t', "timeout", &p.timeout_s, "timeout in seconds"),
        OPT_BOOLEAN('d', "debounds", &p.de_bounds, "use the DE bounds"),
        OPT_BOOLEAN('b', "boolebound", &p.boole_bound, "use the Boole bound"),
        OPT_INTEGER('f', "boolefreq", &p.boole_bound_freq, "freequency at which to use the Boole bound"),
        OPT_END()
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argc = argparse_parse(&argparse, argc, argv);

    if(!p.instance_file) {
        printf("You have to specify an instance file!\n");
        exit(EXIT_FAILURE);
    }

    return p;
}

int main(int argc, const char** argv) {
    Params p = parse_arguments(argc, argv);

    TBKPInstance* instance = tbkp_instance_read(p.instance_file);
    int error = GRBloadenv(&grb_env, NULL);

    if(error) {
        printf("Gurobi loadenv error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    error = GRBsetintparam(grb_env, GRB_INT_PAR_OUTPUTFLAG, 0);

    if(error) {
        printf("Gurobi setintparam output flag error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    TBKPStats stats = tbkp_stats_init(p.timeout_s, p.de_bounds, p.boole_bound, p.boole_bound_freq);
    TBKPSolution* bbsol = tbkp_branch_and_bound(instance, &stats);
    tbkp_stats_to_file(&stats, p.output_file);

    // Clean up
    tbkp_sol_free(&bbsol);
    tbkp_instance_free(&instance); instance = NULL;
    GRBfreeenv(grb_env);

    return 0;
}
