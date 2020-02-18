//
// Created by alberto on 18/02/2020.
//

#include "tbkp_instance.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

TBKPInstance* tbkp_instance_read(const char *const filename) {
    TBKPInstance* instance = malloc(sizeof(*instance));

    if(instance == NULL) {
        printf("Cannot allocate memory for the instance\n");
        exit(EXIT_FAILURE);
    }

    FILE* fd;
    fd = fopen(filename, "r");

    if(fd == NULL) {
        printf("Cannot read %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int n_read;

    n_read = fscanf(fd, "%" SCNuFAST32 " %" SCNuFAST32, &(instance->n_items), &(instance->capacity));

    if(n_read != 2) {
        printf("Error reading the first line with number of items and capacity\n");
        exit(EXIT_FAILURE);
    }

    if(instance->n_items == 0) {
        printf("Error: you have an instance with no items\n");
        exit(EXIT_FAILURE);
    }

    if(instance->n_items > instance->capacity) {
        printf("Error: you have an instance with more items than capacity units\n");
        exit(EXIT_FAILURE);
    }

    instance->weights = malloc(instance->n_items * sizeof(*(instance->weights)));
    instance->profits = malloc(instance->n_items * sizeof(*(instance->profits)));
    instance->probabilities = malloc(instance->n_items * sizeof(*(instance->probabilities)));

    if((instance->weights == NULL) || (instance->profits == NULL) || (instance->probabilities == NULL)) {
        printf("Cannot allocate memory for data\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < instance->n_items; ++i) {
        n_read = fscanf(fd, "%" SCNuFAST32 " %" SCNuFAST32 " %f",
                &(instance->weights[i]), &(instance->profits[i]), &(instance->probabilities[i]));

        if(n_read != 3) {
            printf("Error reading item %zu\n", i);
            exit(EXIT_FAILURE);
        }

        if(instance->weights[i] > instance->capacity) {
            printf("Wrong weight for item %zu: %" PRIuFAST32 " larger than capacity %" PRIuFAST32 "\n",
                    i, instance->weights[i], instance->capacity);
            exit(EXIT_FAILURE);
        }

        if(instance->probabilities[i] < 0.0f || instance->probabilities[i] > 1.0f) {
            printf("Wrong probability for item %zu: %f\n", i, instance->probabilities[i]);
            exit(EXIT_FAILURE);
        }
    }

    if(fclose(fd) != 0) {
        printf("Error closing file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    return instance;
}

void tbkp_instance_free(TBKPInstance** instance_ptr) {
    free((*instance_ptr)->weights); (*instance_ptr)->weights = NULL;
    free((*instance_ptr)->profits); (*instance_ptr)->profits = NULL;
    free((*instance_ptr)->probabilities); (*instance_ptr)->probabilities = NULL;
    free(*instance_ptr); *instance_ptr = NULL;
}

void tbkp_instance_print(const TBKPInstance *const instance) {
    printf("Num items: %" PRIuFAST32 ", Capacity: %" PRIuFAST32 "\n", instance->n_items, instance->capacity);

    for(size_t i = 0; i < instance->n_items; ++i) {
        printf("Object %zu: weight=%" PRIuFAST32 ", profit=%" PRIuFAST32 ", prob=%f\n",
                i, instance->weights[i], instance->profits[i], instance->probabilities[i]);
    }
}