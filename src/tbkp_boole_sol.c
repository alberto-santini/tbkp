//
// Created by alberto on 23/02/2020.
//

#include "tbkp_boole_sol.h"
#include <stdio.h>
#include <stdlib.h>
#include <gurobi_c.h>
#include <assert.h>
#include <time.h>

#define EPS 0.01

TBKPBooleSol tbkp_boolesol_lin_gurobi_get(
        const TBKPInstance* instance,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity)
{
    clock_t start_time = clock();

    GRBmodel* grb_model = NULL;
    int error = GRBnewmodel(grb_env, &grb_model, "booleip", 0, NULL, NULL, NULL, NULL, NULL);
    int n = (int)n_items;

    if(error) {
        printf("Gurobi newmodel error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < n_items; ++i) {
        // x variables:
        error = GRBaddvar(grb_model, 0, NULL, NULL, (double) instance->profits[items[i]], 0.0, 1.0, GRB_BINARY, NULL);

        if (error) {
            printf("Gurobi error adding variable x[%zu]: %d\n", i, error);
            exit(EXIT_FAILURE);
        }
    }

    for(size_t i = 0u; i < n_items; ++i) {
        for(size_t j = 0u; j < n_items; ++j) {
            double coeff = (double) instance->profits[items[i]];
            coeff *= (1.0 - instance->probabilities[items[j]]);
            coeff *= -2.0;  // Double the coefficient because in the model we have the full matrix, not just the
                            // upper triangular one which gurobi requires. So effectively we have both terms
                            // - (1-\pi_{j'}) * p_j  * z_{j j'} [...]
                            // - (1-\pi_j) * p_{j'} * z_{j' j} [...]
                            
            // z variables:
            error = GRBaddvar(grb_model, 0, NULL, NULL, coeff, 0.0, 1.0, GRB_BINARY, NULL);

            if(error) {
                printf("Gurobi error adding variable z[%zu][%zu]: %d\n", i, j, error);
            }
        }
    }

    // 1 constraint: \sum_j w_j x_j <= c
    int* cst_ind = malloc(n_items * sizeof(*cst_ind));
    double* cst_val = malloc(n_items * sizeof(*cst_val));

    if(!cst_ind || !cst_val) {
        printf("Cannot allocate memory for constraint coefficients\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0u; i < n_items; ++i) {
        cst_ind[i] = (int) i;
        cst_val[i] = (double) instance->weights[items[i]];
    }

    error = GRBaddconstr(grb_model, n, cst_ind, cst_val, GRB_LESS_EQUAL, (double)capacity, NULL);

    if(error) {
        printf("Gurobi addconstr error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    // n^2 constraints: z_{ij} >= x_i + x_j - 1
    // i.e. -x_i -x_j + z_{ij} >= -1
    for(size_t i = 0u; i < n_items; ++i) {
        for(size_t j = 0u; j < n_items; ++j) {
            if(i == j) {
                // When i == j, Gurobi doesn't like duplicate indices,
                // so we rewrite the constraint as:
                // z_{ii} = x_i
                // i.e., -x_i + z_{ii} = 0
                int lcst_ind[2] = {
                        (int) i, // x_i
                        n + ((int) i * n + (int) i) // z_{ii}
                };

                double lcst_val[2] = {-1.0, 1.0};

                error = GRBaddconstr(grb_model, 2, lcst_ind, lcst_val, GRB_EQUAL, 0.0, NULL);

                if(error) {
                    printf("Gurobi addcstr error on (%zu, %zu): %d\n", i, j, error);
                    exit(EXIT_FAILURE);
                }
            } else {
                int lcst_ind[3] = {
                        (int) i, // x_i
                        (int) j, // x_j
                        n + ((int)i * n + (int)j) // z_{ij}
                };

                double lcst_val[3] = {-1.0, -1.0, 1.0};

                error = GRBaddconstr(grb_model, 3, lcst_ind, lcst_val, GRB_GREATER_EQUAL, -1.0, NULL);

                if(error) {
                    printf("Gurobi addcstr error on (%zu, %zu): %d\n", i, j, error);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    error = GRBsetintattr(grb_model, "ModelSense", -1);

    if(error) {
        printf("Gurobi setintattr ModelSense error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    error = GRBoptimize(grb_model);

    if(error) {
        printf("Gurobi optimize error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    int grb_status;

    error = GRBgetintattr(grb_model, GRB_INT_ATTR_STATUS, &grb_status);

    if(error) {
        printf("Gurobi error retrieving status: %d\n", error);
        exit(EXIT_FAILURE);
    }

    if(grb_status != GRB_OPTIMAL) {
        printf("Warning: Boole IP status not optimal. It is %d.\n", grb_status);
    }

    double obj;

    error = GRBgetdblattr(grb_model, GRB_DBL_ATTR_OBJVAL, &obj);

    if(error) {
        printf("Gurobi error retrieving obj value: %d\n", error);
        exit(EXIT_FAILURE);
    }

    double* solution = malloc((n_items + n_items * n_items) * sizeof(*solution));

    if(!solution) {
        printf("Cannot allocate memory to retrieve Gurobi solution\n");
        exit(EXIT_FAILURE);
    }

    error = GRBgetdblattrarray(grb_model, GRB_DBL_ATTR_X, 0, n, solution);

    if(error) {
        printf("Gurobi error retrieving solution: %d\n", error);
        exit(EXIT_FAILURE);
    }

    size_t grb_n_items = 0u;

    for(size_t i = 0u; i < n_items; ++i) {
        if(solution[i] > .5) {
            ++grb_n_items;
        }
    }

    size_t* grb_packed_items = malloc(grb_n_items * sizeof(*grb_packed_items));

    if(!grb_packed_items) {
        printf("Cannot allocate memory to save packed objects\n");
        exit(EXIT_FAILURE);
    }

    size_t curr_id = 0u;
    for(size_t i = 0u; i < n_items; ++i) {
        if(solution[i] > .5) {
            grb_packed_items[curr_id++] = items[i];
        }
    }

    GRBfreemodel(grb_model);

    free(cst_ind); cst_ind = NULL;
    free(cst_val); cst_val = NULL;
    free(solution); solution = NULL;

    TBKPBooleSol sol = { .n_items = grb_n_items, .items = grb_packed_items };
    tbkp_boolesol_compute_exact_obj(instance, &sol);

    if(sol.lb < obj - EPS) {
        printf("Error: recomputed objective lower than Boole objective (%f vs %f)\n", sol.lb, obj);
        exit(EXIT_FAILURE);
    }

    clock_t end_time = clock();
    sol.time_to_compute = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    return sol;
}

TBKPBooleSol tbkp_boolesol_quad_gurobi_get(
        const TBKPInstance* instance,
        size_t n_items,
        const size_t* items,
        uint_fast32_t capacity)
{
    clock_t start_time = clock();

    GRBmodel* grb_model = NULL;
    int error = GRBnewmodel(grb_env, &grb_model, "booleqp", 0, NULL, NULL, NULL, NULL, NULL);
    int n = (int)n_items;

    if(error) {
        printf("Gurobi newmodel error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0u; i < n_items; ++i) {
        error = GRBaddvar(grb_model, 0, NULL, NULL, (double) instance->profits[items[i]], 0.0, 1.0, GRB_BINARY, NULL);

        if(error) {
            printf("Gurobi error adding variable %zu: %d\n", i, error);
        }
    }

    /* Quadratic terms:
     *
     * - ( \sum_{i=1}^n p_i x_i ) * ( \sum_{j=1}^n (1-\pi_j) x_j ) =
     *
     * - \sum_{i=1}^n \sum_{j=1}^n p_i (1 - \pi_j) x_i x_j =
     *
     * - \sum_{i=1}^n p_i (1 - \pi_i) x_i^2 -    // Diagonal
     * - \sum_{i=1}^n \sum_{j=i+1}^n [
     *      p_i (1 - \pi_j) + p_j (1 - \pi_i)    // Upper triangular
     *   ] x_i x_j
     */

    int quad_nz = n * (n + 1) / 2;
    int* qrow = malloc((size_t) quad_nz * sizeof(*qrow));
    int* qcol = malloc((size_t) quad_nz * sizeof(*qcol));
    double* qcoeff = malloc((size_t) quad_nz * sizeof(*qcoeff));

    if(!qrow || !qcol || !qcoeff) {
        printf("Cannot allocate memory for quadratic coefficients\n");
        exit(EXIT_FAILURE);
    }

    size_t index = 0u;
    for(size_t i = 0u; i < n_items; ++i) {
        for(size_t j = i; j < n_items; ++j) {
            qrow[index] = (int) i;
            qcol[index] = (int) j;

            if(i == j) {
                qcoeff[index] = - (double)instance->profits[items[i]] * (1 - instance->probabilities[items[i]]);
            } else {
                qcoeff[index] = - (double)instance->profits[items[i]] * (1 - instance->probabilities[items[j]])
                                - (double)instance->profits[items[j]] * (1 - instance->probabilities[items[i]]);
            }

            ++index;
        }
    }

    assert(index == n_items * (n_items + 1) / 2);

    // Quadratic objective coefficients:
    error = GRBaddqpterms(grb_model, quad_nz, qrow, qcol, qcoeff);

    if(error) {
        printf("Gurobi addqpterms error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    int* cst_ind = malloc(n_items * sizeof(*cst_ind));
    double* cst_val = malloc(n_items * sizeof(*cst_val));

    if(!cst_ind || !cst_val) {
        printf("Cannot allocate memory for constraint coefficients\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0u; i < n_items; ++i) {
        cst_ind[i] = (int) i;
        cst_val[i] = (double) instance->weights[items[i]];
    }

    error = GRBaddconstr(grb_model, n, cst_ind, cst_val, GRB_LESS_EQUAL, (double)capacity, NULL);

    if(error) {
        printf("Gurobi addconstr error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    error = GRBsetintattr(grb_model, "ModelSense", -1);

    if(error) {
        printf("Gurobi setintattr ModelSense error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    error = GRBoptimize(grb_model);

    if(error) {
        printf("Gurobi optimize error: %d\n", error);
        exit(EXIT_FAILURE);
    }

    int grb_status;

    error = GRBgetintattr(grb_model, GRB_INT_ATTR_STATUS, &grb_status);

    if(error) {
        printf("Gurobi error retrieving status: %d\n", error);
        exit(EXIT_FAILURE);
    }

    if(grb_status != GRB_OPTIMAL) {
        printf("Warning: Boole BQ status not optimal. It is %d.\n", grb_status);
    }

    double obj;

    error = GRBgetdblattr(grb_model, GRB_DBL_ATTR_OBJVAL, &obj);

    if(error) {
        printf("Gurobi error retrieving obj value: %d\n", error);
        exit(EXIT_FAILURE);
    }

    double* solution = malloc(n_items * sizeof(*solution));

    if(!solution) {
        printf("Cannot allocate memory to retrieve Gurobi solution\n");
        exit(EXIT_FAILURE);
    }

    error = GRBgetdblattrarray(grb_model, GRB_DBL_ATTR_X, 0, n, solution);

    if(error) {
        printf("Gurobi error retrieving solution: %d\n", error);
        exit(EXIT_FAILURE);
    }

    size_t grb_n_items = 0u;

    for(size_t i = 0u; i < n_items; ++i) {
        if(solution[i] > .5) {
            ++grb_n_items;
        }
    }

    size_t* grb_packed_items = malloc(grb_n_items * sizeof(*grb_packed_items));

    if(!grb_packed_items) {
        printf("Cannot allocate memory to save packed objects\n");
        exit(EXIT_FAILURE);
    }

    size_t curr_id = 0u;
    for(size_t i = 0u; i < n_items; ++i) {
        if(solution[i] > .5) {
            grb_packed_items[curr_id++] = items[i];
        }
    }

    GRBfreemodel(grb_model);

    free(qrow); qrow = NULL;
    free(qcol); qcol = NULL;
    free(qcoeff); qcoeff = NULL;
    free(cst_ind); cst_ind = NULL;
    free(cst_val); cst_val = NULL;
    free(solution); solution = NULL;

    TBKPBooleSol sol = { .n_items = grb_n_items, .items = grb_packed_items };
    tbkp_boolesol_compute_exact_obj(instance, &sol);

    if(sol.lb < obj - EPS) {
        printf("Error: recomputed objective lower than Boole objective (%f vs %f)\n", sol.lb, obj);
        exit(EXIT_FAILURE);
    }

    clock_t end_time = clock();
    sol.time_to_compute = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    return sol;
}

void tbkp_boolesol_compute_exact_obj(
        const TBKPInstance *instance,
        TBKPBooleSol *sol)
{
    uint_fast32_t lb_sum = 0.0f;

    for(size_t i = 0u; i < sol->n_items; ++i) {
        lb_sum += instance->profits[sol->items[i]];
    }

    float lb_prod = 1.0f;

    for(size_t i = 0u; i < sol->n_items; ++i) {
        lb_prod *= instance->probabilities[sol->items[i]];
    }

    sol->lb = (float)lb_sum * lb_prod;
    sol->lb_sum_profits = lb_sum;
    sol->lb_product_probabilities = lb_prod;
}

void tbkp_boolesol_print(const TBKPBooleSol *sol) {
    printf("LB value: %f\n", sol->lb);
    printf("Objects packed (%zu):\n\t", sol->n_items);
    for(size_t i = 0; i < sol->n_items; ++i) {
        printf("%zu ", sol->items[i]);
    }
    printf("\n");
}

void tbkp_boolesol_free_inside(TBKPBooleSol* boolesol_ptr) {
    free(boolesol_ptr->items); boolesol_ptr->items = NULL;
}
