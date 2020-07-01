#include "tbkp_dp.h"
#include <combo.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define EPS 1e-6f

static size_t tb_table_index(
    const TBKPInstance *const inst,
    uint_fast32_t max_weight,
    uint_fast32_t min_profit,
    uint_fast32_t max_item)
{
    const uint_fast32_t U = inst->tb_ub_packed_profit;
    return max_weight * (U + 1u) * inst->n_tb_items +
           min_profit * inst->n_tb_items +
           max_item;
}

static float tb_table_get(
    const TBKPInstance *const inst,
    const TBKPDPTables *const t,
    uint_fast32_t max_weight,
    uint_fast32_t min_profit,
    uint_fast32_t max_item)
{
    return t->tb_table[tb_table_index(inst, max_weight, min_profit, max_item)];
}

static void tb_table_set(
    float value,
    const TBKPInstance *const inst,
    TBKPDPTables *const t,
    uint_fast32_t max_weight,
    uint_fast32_t min_profit,
    uint_fast32_t max_item)
{
    t->tb_table[tb_table_index(inst, max_weight, min_profit, max_item)] = value;
}

static size_t det_table_index(
    const TBKPInstance *const inst,
    uint_fast32_t max_weight,
    uint_fast32_t max_item)
{
    // In the instance deterministic items are indexed from
    // inst->last_tb_item +1u to inst->n_items - 1u. In the
    // DP table, however, we index the deterministic items
    // from 0 to inst->n_det_items - 1u, so we have to convert
    // the indices.
    const uint_fast32_t det_item = inst->last_tb_item_index + max_item + 1u;

    return max_weight * inst->n_det_items +
           det_item;
}

static uint_fast32_t det_table_get(
    const TBKPInstance *const inst,
    const TBKPDPTables *const t,
    uint_fast32_t max_weight,
    uint_fast32_t max_item)
{
    return t->det_table[det_table_index(inst, max_weight, max_item)];
}

static void det_table_set(
    uint_fast32_t value,
    const TBKPInstance *const inst,
    TBKPDPTables *const t,
    uint_fast32_t max_weight,
    uint_fast32_t max_item)
{
    t->det_table[det_table_index(inst, max_weight, max_item)] = value;
}

static float tb_solve_recursive(
    uint_fast32_t d, uint_fast32_t v, uint_fast32_t j,
    const TBKPInstance *const inst, TBKPDPTables* t)
{
    const float val = tb_table_get(inst, t, d, v, j);
    if(val != TB_TABLE_UNKNOWN) {
        return val;
    }

    assert(j > 0u);

    const float first_val = tb_solve_recursive(d, v, j - 1u, inst, t);

    float second_val;
    if(d < inst->weights[j] || v < inst->profits[j]) {
        second_val = 0.0f;
    } else {
        second_val = tb_solve_recursive(
            d - inst->weights[j],
            v - inst->profits[j],
            j - 1u, inst, t) * inst->probabilities[j];
    }

    const float newval = (first_val > second_val) ? first_val : second_val;
    tb_table_set(newval, inst, t, d, v, j);

    return newval;
}

static uint_fast32_t det_solve_recursive(
    uint_fast32_t d, uint_fast32_t j,
    const TBKPInstance* inst, TBKPDPTables* t)
{
    const uint_fast32_t val = det_table_get(inst, t, d, j);
    if(val != DET_TABLE_UNKNOWN) {
        return val;
    }

    assert(j > 0u);

    const uint_fast32_t first_val = det_solve_recursive(d, j - 1u, inst, t);

    uint_fast32_t second_val;
    const uint_fast32_t wj = inst->weights[inst->last_tb_item_index + j + 1u];
    if(d < wj) {
        second_val = 0.0f;
    } else {
        second_val = det_solve_recursive(d - wj, j - 1u, inst, t);
    }

    const uint_fast32_t newval = (first_val > second_val) ? first_val : second_val;
    det_table_set(newval, inst, t, d, j);

    return newval;
}

static uint_fast32_t det_solve_fixed(
    uint_fast32_t d,
    const TBKPInstance *const inst,
    TBKPDPTables* t)
{
    return det_solve_recursive(d, inst->n_det_items, inst, t);
}

static float tb_solve_fixed(
    uint_fast32_t d,
    uint_fast32_t v,
    const TBKPInstance *const inst,
    TBKPDPTables* t)
{
    return tb_solve_recursive(d, v, inst->n_tb_items, inst, t);
}

float tbkp_dp_solve(const TBKPInstance *const inst) {
    const uint_fast32_t U = inst->tb_ub_packed_profit;
    TBKPDPTables t = tbkp_dp_tables_init(inst);
    float best_sol = -1.0f;

    for(uint_fast32_t d = 0; d <= inst->capacity; ++d) {
        for(uint_fast32_t v = 0; v <= U; ++v) {
            const float tb_sol = tb_solve_fixed(d, v, inst, &t);
            const uint_fast32_t det_sol = det_solve_fixed(inst->capacity - d, inst, &t);
            const float sol = (float)(v + det_sol) * tb_sol;

            if(sol > best_sol) {
                best_sol = sol;
            }
        }
    }

    return best_sol;
}

TBKPDPTables tbkp_dp_tables_init(const TBKPInstance *const inst) {
    TBKPDPTables t;

    const size_t tb_table_sz =
        (inst->capacity + 1u) *
        (inst->tb_ub_packed_profit + 1u) *
        inst->n_tb_items;

    t.tb_table = calloc(tb_table_sz, sizeof(*t.tb_table));

    if(!t.tb_table) {
        printf("Not enough memory for the TB items DP table!\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0u; i < tb_table_sz; ++i) {
        t.tb_table[i] = TB_TABLE_UNKNOWN;
    }

    for(uint_fast32_t d = 0u; d <= inst->capacity; ++d) {
        tb_table_set(1.0f, inst, &t, d, 0u, 0u);

        for(uint_fast32_t v = 0u; v <= inst->tb_ub_packed_profit; ++v) {
            tb_table_set(0.0f, inst, &t, d, v, 0u);
        }
    }

    const size_t det_table_sz = (inst->capacity + 1u) * inst->n_det_items;

    t.det_table = calloc(det_table_sz, sizeof(*t.det_table));

    if(!t.det_table) {
        printf("Not enough memory for the DET items DP table!\n");
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0u; i < det_table_sz; ++i) {
        t.det_table[i] = DET_TABLE_UNKNOWN;
    }

    for(uint_fast32_t d = 0u; d <= inst->capacity; ++d) {
        det_table_set(0.0f, inst, &t, d, 0u);
    }

    for(uint_fast32_t j = 0u; j < inst->n_det_items; ++j) {
        for(uint_fast32_t d = 0u; d < inst->weights[inst->last_tb_item_index + j + 1u]; ++d) {
            det_table_set(0.0f, inst, &t, d, j);
        }
        for(uint_fast32_t d = inst->weights[inst->last_tb_item_index + j + 1u]; d <= inst->capacity; ++d) {
            det_table_set(inst->profits[inst->last_tb_item_index + j + 1u], inst, &t, d, j);
        }
    }

    return t;
}

void tbkp_dp_tables_free(TBKPDPTables* t) {
    free(t->tb_table); t->tb_table = NULL;
    free(t->det_table); t->det_table = NULL;
}