#include "tbkp_cr_sol.h"
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

#define EPSC 1e-6

/************************************************************
 * LOCAL HELPER FUNCTIONS                                   *
 ************************************************************/

static float bin_search(size_t n, float* x, float* d, uint_fast32_t *p, float* q) {
    // Next point is z = x + lambda *x   with lambda \in [0,1]
    float bestobj = 0.0f;
    float bestl = -1.0f;

    float ptot = 0.0f;
    float product = 1.0f;

    // Iterative loop
    float l0 = 0.0f;
    float l1 = 1.0f;
    while(l1 - l0 > EPSC) {
        // Define two lambda values in the current interval
        float lgap = l1 - l0;
        float delta = lgap / 3;
        float beta = l0 + delta;
        float gamma = l0 + 2 * delta;

        // Evaluate the first point
        ptot = 0.0f;
        product = 1.0f;
        for(size_t j = 0; j < n; j++) {
            const float zj = x[j] + beta * d[j];
            ptot += (float)(p[j]) * zj;
            const float alpha = 1.0f - q[j] * zj;
            product = product * alpha;
        }

        float objbeta = ptot * product;
        if(objbeta > bestobj) {
                bestobj = objbeta;
                bestl = beta;
        }

        // Evaluate the second point
        ptot = 0.0f;
        product = 1.0f;
        for(size_t j = 0; j < n; j++) {
            const float zj = x[j] + gamma * d[j];
            ptot += (float)(p[j]) * zj;
            const float alpha = 1.0f - q[j] * zj;
            product = product * alpha;
        }

        float objgamma = ptot * product;
        if(objgamma > bestobj) {
            bestobj = objgamma;
            bestl = gamma;
        }

        if(objbeta < objgamma) {
            // remove interval [l0, beta]
            l0 = beta;
        } else {
            // remove interval [gamma, l1]
            l1 = gamma;
        }
    }

    return bestl;
}

static void solve_Dantzig(size_t n, float* p, uint_fast32_t* w, uint_fast32_t c, const TBKPBBFixedStatus* xbra, float* y) {
    uint_fast32_t cres = c;
    int* flag = (int*) calloc(n, sizeof(int));

    for(size_t j = 0; j < n; j++) {
        y[j] = 0;
        flag[j] = 0;

        if(xbra[j] != -1) {
            flag[j] = 1;
            
            if(xbra[j] == 1) {
                y[j] = 1.0;
                cres -= w[j];
            }
        }
    }

    for(size_t cont = 0; cont < n; cont++) {
        if(cres <= 0) {
            break;
        }

        float maxratio = 0.0;
        int index = -1;
        for(size_t j = 0; j < n; j++) {
            if(flag[j] == 0) {
                const float ratio = (float)(p[j]) / (float)(w[j]);
                if(ratio > maxratio) {
                    maxratio = ratio;
                    index = (int)j;
                }
            }
        }

        if(index < 0) {
            break;
        }

        flag[index] = 1;
        if(cres > w[index]) {
            y[index] = 1.0f;
            cres -= w[index];
        } else {
            y[index] = 1.0f * (float)(cres) / (float)(w[index]);
            cres = 0;
        }
    }

    free(flag);
}

TBKPContinuousRelaxationSol tbkp_crsol_get(
        const TBKPInstance *const instance,
        const TBKPBBFixedStatus *const xbra)
{
    clock_t start_time = clock();

    size_t n = instance->n_items;
    float* xc = (float*) calloc(n, sizeof(float)); // x variables
    float* a = (float*) calloc(n, sizeof(float)); // \alpha variables
    float* q = (float*) calloc(n, sizeof(float)); // q probabilities

    for(size_t j = 0; j < n; j++) {
        if(xbra[j] == 1) {
            xc[j] = 1.0;
        } else {
            xc[j] = 0.0;
        }

        a[j] = 1.0f;
        q[j] = 1.0f - instance->probabilities[j];
    }

    float product = 1.0;
    float ptot = 0.0;
    float* g = (float*) calloc(n, sizeof(float)); // gradient
    float* d = (float*) calloc(n, sizeof(float)); // improving direction (if any)
    float* y = (float*) calloc(n, sizeof(float)); // next solution
    float objval = 0.0;

    while(true) {
        float oldobj = objval;

        // Compute the gradient
        product = 1.0;
        for(size_t j = 0; j < n; j++) {
            product *= (1.0f - q[j]*xc[j]);
        }
        for(size_t j = 0; j < n; j++) {
            const float value = (float)(instance->profits[j]) - ptot * q[j] / (1.0f - q[j] * xc[j]);
            g[j] = product * value;
        }
        
        // Find a candidate direction for improving
        solve_Dantzig(n, g, instance->weights, instance->capacity, xbra, y);
        for(size_t j = 0; j < n; j++) {
            d[j] = y[j] - xc[j];
        }

        // Check if the direction is improving
        float delta = 0.0f;
        for(size_t j = 0; j < n; j++) {
            delta += g[j] * d[j];
        }
        if(delta < EPSC) {
            break;
        }

        // Determine the optimal value for parameter \lambda
        float lambda = bin_search(n, xc, d, instance->profits, q);

        // Determine the next point
        ptot = 0.0f;
        product = 1.0f;
        for(size_t j = 0; j < n; j++) {
            const float zj = xc[j] +  lambda * d[j];
            ptot += (float)(instance->profits[j]) * zj;
            a[j] = 1.0f - q[j] * zj;
            product = product * a[j];
            xc[j] = zj;
        }
        objval = ptot * product;
        
        // Break if no improvement
        if(objval - oldobj < EPSC) {
            break;
        }
    }

    _Bool integer_sol = true;
    size_t n_items = 0u;
    for(size_t i = 0u; i < n; ++i) {    
        if((xc[i] > EPSC) && (xc[i] < 1.0f - EPSC)) {
            integer_sol = false;
        }

        if(xc[i] > 1.0f - EPSC) {
            ++n_items;
        }
    }

    size_t* items = NULL;

    if(integer_sol) {
        items = (size_t*) calloc(n_items, sizeof(*items));

        size_t id = 0u;
        for(size_t i = 0u; i < n; ++i) {
            if(xc[i] > 1.0f - EPSC) {
                assert(id < n_items);
                items[id++] = i;
            }
        }
    }

    free(y);
    free(d);
    free(g);
    free(q);
    free(a);
    free(xc);

    clock_t end_time = clock();
    float elapsed_time = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    return (TBKPContinuousRelaxationSol){
        .ub = objval,
        .is_lb = integer_sol,
        .lb_product_probabilities = product,
        .lb_sum_profits = (uint_fast32_t)ptot,
        .n_items = n_items,
        .items = items,
        .time_to_compute = elapsed_time
    };
}

void tbkp_crsol_free_inside(TBKPContinuousRelaxationSol *const crsol_ptr) {
    if(crsol_ptr->items) {
        free(crsol_ptr->items); crsol_ptr->items = NULL;
    }
}

void tbkp_crsol_print(const TBKPContinuousRelaxationSol *const sol) {
    printf("UB value: %f\n", sol->ub);
    printf("Solution integer?: %d\n", sol->is_lb);
}
