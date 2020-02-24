#include "tbkp_instance.h"
#include "tbkp_de_sol.h"
#include "tbkp_boole_sol.h"
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

    tbkp_instance_print(instance);

    size_t* items = malloc(instance->n_items * sizeof(*items));
    for(size_t i = 0; i < instance->n_items; ++i) {
        items[i] = i;
    }

    TBKPDeterministicEqSol desol = tbkp_desol_get(instance, instance->n_items, items, instance->capacity);

    tbkp_desol_print(&desol);
    tbkp_desol_free_inside(&desol);

    TBKPBoolSol boolsol = tbkp_boolsol_quad_gurobi_get(instance, instance->n_items, items, instance->capacity);

    tbkp_boolsol_print(&boolsol);
    printf("Checking if exact objective value is better...\n");
    tbkp_boolsol_compute_exact_obj(instance, &boolsol);
    tbkp_boolsol_print(&boolsol);
    tbkp_boolsol_free_inside(&boolsol);

    GRBfreeenv(grb_env);

    tbkp_instance_free(&instance);

    free(items); items = NULL;

    return 0;
}
