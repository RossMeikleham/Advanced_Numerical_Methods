#ifndef P_PIVOT_GAUSSIAN_ELIM_H
#define P_PIVOT_GAUSSIAN_ELIM_H

#include "../matrix.h"
#include <iostream>
#include <vector>
#include <tuple>


template <class Type>
matrix<Type> pivoting_gauss_elimination(matrix<Type> a1) {
    auto a = a1;
    for (uint64_t i = 0; i < a.n() - 1; i++) {
        
        auto r = i;
        auto max = abs(a.get(i, i));
        for (uint64_t j = i + 1; j < a.m(); j++) {
            if (abs(a.get(j, i)) > max) {
                r = j;
                max = abs(a.get(j, i));
            }
        }
        
        auto temp = a.row(r);
        a.setRow(r, a.row(i));
        a.setRow(i, temp);
        
        #pragma omp parallel for
        for (uint64_t j = i + 1; j < a.m(); j++) {
            a.setRow(j, a.row(j) - ((a.get(j, i) * a.row(i)) / a.get(i, i)));          
        }
    }

    return a;
}

#endif
