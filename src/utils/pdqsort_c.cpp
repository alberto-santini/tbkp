//
// Created by alberto on 20/02/2020.
//

#include "pdqsort.h"
#include "combo.h"
#include <cstddef>

extern "C" {
    void pdqsort_c(item* ary, size_t n) {
        pdqsort(ary, ary + n, [](const auto &it1, const auto &it2) -> bool { return it1.id < it2.id; });
    }
}