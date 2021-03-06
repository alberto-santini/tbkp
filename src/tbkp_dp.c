#include "tbkp_dp.h"
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

static void tbkp_dp_timebombtable_set_on_active_row(
    uint_fast32_t d, uint_fast32_t v, double val,
    TBKPDPTimeBombTable *const tb_t,
    const TBKPInstance *const inst
) {
    const uint_fast32_t U = inst->tb_ub_packed_profit;
    const size_t index = d * (U + 1u) + v;

    assert(index < (inst->capacity + 1u) * (U + 1u));
    assert(tb_t->active_row <= 1);

    tb_t->rows[tb_t->active_row][index] = val;
}

static void tbkp_dp_timebombtable_set_on_new_row(
    uint_fast32_t d, uint_fast32_t v, double val,
    TBKPDPTimeBombTable *const tb_t,
    const TBKPInstance *const inst
) {
    const uint_fast32_t U = inst->tb_ub_packed_profit;
    const size_t index = d * (U + 1u) + v;

    assert(index < (inst->capacity + 1u) * (U + 1u));
    assert(tb_t->active_row <= 1);

    tb_t->rows[1u - tb_t->active_row][index] = val;
}

static double tbkp_dp_timebombtable_get(
    uint_fast32_t d, uint_fast32_t v, uint_fast32_t j,
    const TBKPDPTimeBombTable *const tb_t,
    const TBKPInstance *const inst
) {
    assert(j == tb_t->last_index || j == tb_t->last_index - 1u);

    const size_t row_idx = (j == tb_t->last_index) ? tb_t->active_row : (1u - tb_t->active_row);
    const uint_fast32_t U = inst->tb_ub_packed_profit;
    const size_t index = d * (U + 1u) + v;

    return tb_t->rows[row_idx][index];
}

static double tbkp_dp_timebombtable_getopt(
    uint_fast32_t d, uint_fast32_t v,
    const TBKPDPTimeBombTable *const tb_t,
    const TBKPInstance *const inst
) {
    return tbkp_dp_timebombtable_get(d, v, inst->n_tb_items - 1u, tb_t, inst);
}

// Returns true if it times out - otherwise, false.
static _Bool tbkp_dp_timebombtable_compute(
    TBKPDPTimeBombTable* tb_t,
    clock_t start_time,
    const TBKPInstance *const inst,
    const TBKPParams *const params
) {
    assert(tb_t->last_index == 0u);
    const uint_fast32_t c = inst->capacity;
    const uint_fast32_t U = inst->tb_ub_packed_profit;

    for(uint_fast32_t j = 1u; j <= inst->n_tb_items; ++j) {
        const uint_fast32_t wj = inst->weights[j];
        const uint_fast32_t pj = inst->profits[j];
        const double pij = inst->probabilities[j];

        for(uint_fast32_t dd = 0u; dd <= c; ++dd) {
            const uint_fast32_t d = c - dd;

            for(uint_fast32_t v = 0u; v <= U; ++v) {
                const double current_val = tbkp_dp_timebombtable_get(d, v, j - 1u, tb_t, inst);

                if(wj > d || pj > v) {
                    tbkp_dp_timebombtable_set_on_new_row(d, v, current_val, tb_t, inst);
                } else {
                    const double new_val = tbkp_dp_timebombtable_get(d - wj, v - pj, j - 1u, tb_t, inst) * pij;

                    if(new_val > current_val) {
                        tbkp_dp_timebombtable_set_on_new_row(d, v, new_val, tb_t, inst);
                    } else {
                        tbkp_dp_timebombtable_set_on_new_row(d, v, current_val, tb_t, inst);
                    }
                }

                const clock_t current_time = clock();

                if((float)(current_time - start_time) / CLOCKS_PER_SEC > params->timeout) {
                    return true;
                }
            }
        }

        if(j < inst->n_tb_items) {
            ++(tb_t->last_index);
            tb_t->active_row = 1u - tb_t->active_row;
        }
    }

    return false;
}

static void tbkp_dp_deterministictable_set_raw(
    uint_fast32_t d, uint_fast32_t j, uint_fast32_t val,
    TBKPDPDeterministicTable* d_t,
    const TBKPInstance *const inst
) {
    const uint_fast32_t ndet = inst->n_det_items;
    const size_t index = d * (ndet + 1u) + j;

    assert(index < (inst->capacity + 1u) * (ndet + 1u));
    
    d_t->t[index] = val;
}

static void tbkp_dp_deterministictable_set(
    uint_fast32_t d, uint_fast32_t j, uint_fast32_t val,
    TBKPDPDeterministicTable* d_t,
    const TBKPInstance *const inst
) {
    tbkp_dp_deterministictable_set_raw(d, j + 1u, val, d_t, inst);
}

static uint_fast32_t tbkp_dp_deterministictable_get(
    uint_fast32_t d, uint_fast32_t j,
    const TBKPDPDeterministicTable *const d_t,
    const TBKPInstance *const inst
) {
    j = j + 1u;
    const uint_fast32_t ndet = inst->n_det_items;
    const size_t index = d * (ndet + 1u) + j;

    assert(index < (inst->capacity + 1u) * (ndet + 1u));

    return d_t->t[index];
}

static uint_fast32_t tbkp_dp_deterministictable_getopt(
    uint_fast32_t d,
    const TBKPDPDeterministicTable *const d_t,
    const TBKPInstance *const inst
) {
    return tbkp_dp_deterministictable_get(d, inst->n_det_items - 1u, d_t, inst);
}

static _Bool tbkp_dp_deterministictable_compute(
    TBKPDPDeterministicTable* d_t,
    clock_t start_time,
    const TBKPInstance *const inst,
    const TBKPParams *const params
) {

    for(uint_fast32_t j = 1u; j < inst->n_det_items; ++j) {
        const uint_fast32_t j_idx = inst->last_tb_item_index + j + 1u;
        const uint_fast32_t wj = inst->weights[j_idx];
        const uint_fast32_t pj = inst->profits[j_idx];

        for(uint_fast32_t dd = 0u; dd <= inst->capacity; ++dd) {
            const uint_fast32_t d = inst->capacity - dd;
            const uint_fast32_t old_value = tbkp_dp_deterministictable_get(d, j - 1u, d_t, inst);
            uint_fast32_t new_value = 0;

            if(d >= wj) {
                new_value = pj + tbkp_dp_deterministictable_get(d - wj, j - 1u, d_t, inst);
            }

            if(new_value > old_value) {
                tbkp_dp_deterministictable_set(d, j, new_value, d_t, inst);
            } else {
                tbkp_dp_deterministictable_set(d, j, old_value, d_t, inst);
            }

            const clock_t current_time = clock();

            if((float)(current_time - start_time) / CLOCKS_PER_SEC > params->timeout) {
                return true;
            }
        }
    }

    return false;
}

TBKPDPTimeBombTable tbkp_dp_timebombtable_init(const TBKPInstance* inst) {
    TBKPDPTimeBombTable tb_t;
    const uint_fast32_t c = inst->capacity;
    const uint_fast32_t U = inst->tb_ub_packed_profit;

    tb_t.active_row = 0u;
    tb_t.last_index = 0u;
    tb_t.rows[0] = calloc((c + 1u) * (U + 1u), sizeof(*(tb_t.rows[0])));
    tb_t.rows[1] = calloc((c + 1u) * (U + 1u), sizeof(*(tb_t.rows[1])));

    if(!tb_t.rows[0] || !tb_t.rows[1]) {
        printf("Cannot allocate memory for TB Dynamic Programming table!\n");
        exit(EXIT_FAILURE);
    }

    for(uint_fast32_t d = 0; d <= c; ++d) {
        tbkp_dp_timebombtable_set_on_active_row(d, 0, 1.0, &tb_t, inst);
    }

    for(uint_fast32_t d = inst->weights[0u]; d <= c; ++d) {
        tbkp_dp_timebombtable_set_on_active_row(d, inst->profits[0], inst->probabilities[0], &tb_t, inst);
    }

    return tb_t;
}

void tbkp_dp_timebombtable_free(TBKPDPTimeBombTable* tb_t) {
    free(tb_t->rows[0]); tb_t->rows[0] = NULL;
    free(tb_t->rows[1]); tb_t->rows[1] = NULL;
}

TBKPDPDeterministicTable tbkp_dp_deterministictable_init(const TBKPInstance *const inst) {
    TBKPDPDeterministicTable d_t;
    const uint_fast32_t c = inst->capacity;
    const uint_fast32_t n = inst->n_det_items;

    d_t.t = calloc((c + 1u) * (n + 1u), sizeof(*(d_t.t)));

    for(uint_fast32_t d = 0u; d <= c; ++d) {
        tbkp_dp_deterministictable_set_raw(d, 0u, 0u, &d_t, inst);
    }

    const uint_fast32_t first_det_obj_weight = inst->weights[inst->last_tb_item_index + 1u];
    const uint_fast32_t first_det_obj_profit = inst->profits[inst->last_tb_item_index + 1u];

    for(uint_fast32_t d = 0u; d < first_det_obj_weight; ++d) {
        tbkp_dp_deterministictable_set(d, 0u, 0u, &d_t, inst);
    }

    for(uint_fast32_t d = first_det_obj_weight; d <= c; ++d) {
        tbkp_dp_deterministictable_set(d, 0u, first_det_obj_profit, &d_t, inst);
    }

    return d_t;
}

void tbkp_dp_deterministictable_free(TBKPDPDeterministicTable* d_t) {
    free(d_t->t); d_t->t = NULL;
}

double tbkp_dp_solve(const TBKPInstance *const inst, TBKPDPStats* stats, const TBKPParams *const params) {
    const uint_fast32_t U = inst->tb_ub_packed_profit;
    TBKPDPTimeBombTable tb_t = tbkp_dp_timebombtable_init(inst);
    TBKPDPDeterministicTable d_t = tbkp_dp_deterministictable_init(inst);

    const clock_t start_time = clock();
    _Bool timeout = false;
    stats->elapsed_time = params->timeout;
    stats->built_timebomb_table = false;
    stats->built_deterministic_table = false;

    timeout = tbkp_dp_timebombtable_compute(&tb_t, start_time, inst, params);

    if(timeout) {
        return -1.0;
    }

    stats->built_timebomb_table = true;

    timeout = tbkp_dp_deterministictable_compute(&d_t, start_time, inst, params);

    if(timeout) {
        return -1.0;
    }

    stats->built_deterministic_table = true;

    double best_sol = -1.0;

    for(uint_fast32_t d = 0; d <= inst->capacity; ++d) {
        const uint_fast32_t det_sol = tbkp_dp_deterministictable_getopt(inst->capacity - d, &d_t, inst);

        for(uint_fast32_t v = 0; v <= U; ++v) {
            double sol = ((double)v + (double)det_sol);
            sol *= tbkp_dp_timebombtable_getopt(d, v, &tb_t, inst);

            if(sol > best_sol) {
                best_sol = sol;
            }
        }
    }

    const clock_t end_time = clock();
    stats->elapsed_time = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    tbkp_dp_timebombtable_free(&tb_t);
    tbkp_dp_deterministictable_free(&d_t);

    return best_sol;
}
