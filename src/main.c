#include "tbkp_instance.h"
#include "tbkp_ub.h"
#include <stddef.h>
#include <stdlib.h>

int main() {
    TBKPInstance* instance = tbkp_instance_read("../data/example.txt");

    tbkp_instance_print(instance);

    size_t* items = malloc(instance->n_items * sizeof(*items));
    for(size_t i = 0; i < instance->n_items; ++i) {
        items[i] = i;
    }

    TBKPDeterministicEqUB ub = tbkp_deub_get(instance, instance->n_items, items, instance->capacity);

    tbkp_deub_print(&ub);

    tbkp_deub_free_inside(&ub);
    tbkp_instance_free(&instance);

    return 0;
}
