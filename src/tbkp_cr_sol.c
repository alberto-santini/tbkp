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

static double bin_search(size_t n, double* x, double* d, uint_fast32_t *p, double* q) {
    // Next point is z = x + lambda *x   with lambda \in [0,1]
    double bestobj = 0.0;
    double bestl = -1.0;

    double ptot = 0.0;
    double product = 1.0;

    // Iterative loop
    double l0 = 0.0;
    double l1 = 1.0;
    while(l1 - l0 > EPSC) {
        // Define two lambda values in the current interval
        double lgap = l1 - l0;
        double delta = lgap / 3;
        double beta = l0 + delta;
        double gamma = l0 + 2 * delta;

        // Evaluate the first point
        ptot = 0.0;
        product = 1.0;
        for(size_t j = 0; j < n; j++) {
            const double zj = x[j] + beta * d[j];
            ptot += (double)(p[j]) * zj;
            const double alpha = 1.0 - q[j] * zj;
            product = product * alpha;
        }

        double objbeta = ptot * product;
        if(objbeta > bestobj) {
                bestobj = objbeta;
                bestl = beta;
        }

        // Evaluate the second point
        ptot = 0.0;
        product = 1.0;
        for(size_t j = 0; j < n; j++) {
            const double zj = x[j] + gamma * d[j];
            ptot += (double)(p[j]) * zj;
            const double alpha = 1.0 - q[j] * zj;
            product = product * alpha;
        }

        double objgamma = ptot * product;
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

static void solve_Dantzig(size_t n, double* p, uint_fast32_t* w, uint_fast32_t c, const TBKPBBFixedStatus* xbra, double* y) {
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

        double maxratio = 0.0;
        int index = -1;
        for(size_t j = 0; j < n; j++) {
            if(flag[j] == 0) {
                const double ratio = (double)(p[j]) / (double)(w[j]);
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
            y[index] = 1.0;
            cres -= w[index];
        } else {
            y[index] = 1.0 * (double)(cres) / (double)(w[index]);
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
    double* xc = (double*) calloc(n, sizeof(double)); // x variables
    double* a = (double*) calloc(n, sizeof(double)); // \alpha variables
    double* q = (double*) calloc(n, sizeof(double)); // q probabilities

    for(size_t j = 0; j < n; j++) {
        if(xbra[j] == 1) {
            xc[j] = 1.0;
        } else {
            xc[j] = 0.0;
        }

        a[j] = 1.0;
        q[j] = 1.0 - instance->probabilities[j];
    }

    double product = 1.0;
    double ptot = 0.0;
    double* g = (double*) calloc(n, sizeof(double)); // gradient
    double* d = (double*) calloc(n, sizeof(double)); // improving direction (if any)
    double* y = (double*) calloc(n, sizeof(double)); // next solution
    double objval = 0.0;

    while(true) {
        double oldobj = objval;

        // Compute the gradient
        product = 1.0;
        for(size_t j = 0; j < n; j++) {
            product *= (1.0 - q[j]*xc[j]);
        }
        for(size_t j = 0; j < n; j++) {
            const double value = (double)(instance->profits[j]) - ptot * q[j] / (1.0 - q[j] * xc[j]);
            g[j] = product * value;
        }
        
        // Find a candidate direction for improving
        solve_Dantzig(n, g, instance->weights, instance->capacity, xbra, y);
        for(size_t j = 0; j < n; j++) {
            d[j] = y[j] - xc[j];
        }

        // Check if the direction is improving
        double delta = 0.0;
        for(size_t j = 0; j < n; j++) {
            delta += g[j] * d[j];
        }
        if(delta < EPSC) {
            break;
        }

        // Determine the optimal value for parameter \lambda
        double lambda = bin_search(n, xc, d, instance->profits, q);

        // Determine the next point
        ptot = 0.0;
        product = 1.0;
        for(size_t j = 0; j < n; j++) {
            const double zj = xc[j] +  lambda * d[j];
            ptot += (double)(instance->profits[j]) * zj;
            a[j] = 1.0 - q[j] * zj;
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
        if((xc[i] > EPSC) && (xc[i] < 1.0 - EPSC)) {
            integer_sol = false;
        }

        if(xc[i] > 1.0 - EPSC) {
            ++n_items;
        }
    }

    size_t* items = NULL;

    if(integer_sol) {
        items = (size_t*) calloc(n_items, sizeof(*items));

        size_t id = 0u;
        for(size_t i = 0u; i < n; ++i) {
            if(xc[i] > 1.0 - EPSC) {
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
