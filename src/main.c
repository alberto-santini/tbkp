#include "tbkp_instance.h"
#include "tbkp_de_sol.h"
#include "tbkp_boole_sol.h"
#include "tbkp_BaB.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <gurobi_c.h>

GRBenv* grb_env = NULL;

int main() {
    TBKPInstance* instance = tbkp_instance_read("../data/example.txt");
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

    printf("Instance:\n");
    tbkp_instance_print(instance);
    printf("\n\n");

    size_t* items = malloc(instance->n_items * sizeof(*items));
    for(size_t i = 0; i < instance->n_items; ++i) {
        items[i] = i;
    }

    printf("Getting an upper bound solving a deterministic 01KP\n");
    TBKPDeterministicEqSol desol = tbkp_desol_get(instance, instance->n_items, items, instance->capacity);

    tbkp_desol_print(&desol);
    tbkp_desol_free_inside(&desol);

    printf("Getting LB solving the quadratic programme...\n");
    TBKPBoolSol boolsol = tbkp_boolsol_quad_gurobi_get(instance, instance->n_items, items, instance->capacity);

    tbkp_boolsol_print(&boolsol);
    printf("Checking if exact objective value is better...\n");
    tbkp_boolsol_compute_exact_obj(instance, &boolsol);
    tbkp_boolsol_print(&boolsol);
    tbkp_boolsol_free_inside(&boolsol);

    printf("Getting LB solving the linearisation...\n");
    TBKPBoolSol lin_boolsol = tbkp_boolsol_lin_gurobi_get(instance, instance->n_items, items, instance->capacity);

    tbkp_boolsol_print(&lin_boolsol);
    printf("Checking if the exact objective value is better...\n");
    tbkp_boolsol_compute_exact_obj(instance, &lin_boolsol);
    tbkp_boolsol_print(&lin_boolsol);
    tbkp_boolsol_free_inside(&lin_boolsol);

    printf("Launching the branch-and-bound algorithm...\n");
    TBKPSolution* bbsol = tbkp_branch_and_bound(instance);

    tbkp_sol_print(bbsol, instance);
    tbkp_sol_free(&bbsol);

    // Clean up
    GRBfreeenv(grb_env);
    tbkp_instance_free(&instance);
    free(items); items = NULL;

    return 0;
}
