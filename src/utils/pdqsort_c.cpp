//
// Created by alberto on 20/02/2020.
//

#include "pdqsort.h"
#include <cstddef>

extern "C" {
    void pdqsort_c(size_t* ary, size_t n) {
        pdqsort(ary, ary + n);
    }
}