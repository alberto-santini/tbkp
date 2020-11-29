//
// Created by alberto on 21/04/2020.
//

#include "tbkp_bb_solution.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

TBKPBBSolution* tbkp_sol_init(const TBKPInstance *const instance) {
    TBKPBBSolution* solution = malloc(sizeof(*solution));

    if(!solution) {
        printf("Cannot allocate memory for solution\n");
        exit(EXIT_FAILURE);
    }

    solution->x = malloc(instance->n_items * sizeof(*solution->x));

    if(!solution->x) {
        printf("Cannot allocate memory for solution's x vector\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0u; i < instance->n_items; ++i) {
        solution->x[i] = false;
    }

    solution->prod_probabilities = 1.0;
    solution->sum_profits = 0u;
    solution->value = 0.0;

    return solution;
}

void tbkp_sol_print(const TBKPBBSolution *const solution, const TBKPInstance *const instance) {
    printf("Solution with value %.2f (%" PRIuFAST32 " sum of profits, %.2f prod of probabilities)\n",
            solution->value, solution->sum_profits, solution->prod_probabilities);
    printf("Packed objects:\n");

    for(size_t i = 0u; i < instance->n_items; ++i) {
        if(solution->x[i]) {
            printf("\tObj %zu, profit %3" PRIuFAST32 ", weight %3" PRIuFAST32 ", prob %.2f\n",
                    i, instance->profits[i], instance->weights[i], instance->probabilities[i]);
        }
    }
}

void tbkp_sol_free(TBKPBBSolution** solution_ptr) {
    free((*solution_ptr)->x); (*solution_ptr)->x = NULL;
    free(*solution_ptr); *solution_ptr = NULL;
}