#ifndef GAUSSIAN_ELIM_H
#define GAUSSIAN_ELIM_H

#include "../matrix.h"
#include <iostream>
#include <vector>
#include <tuple>


template <class Type>
matrix<Type> naive_gauss_elimination(matrix<Type> const ax) {
    auto ax2 = ax;
    
    for (uint64_t k = 0; k < ax.m() - 1; k++) {
        for (uint64_t j = k + 1; j < ax.m(); j++) {
            ax2.set(j, k, ax2.get(j, k) / ax2.get(k, k));
            
            #pragma omp parallel for 
            for (uint64_t i = k + 1; i < ax2.m(); i++) {
                ax2.set(j, i, ax2.get(j, i) - (ax2.get(j, k) * ax2.get(k, i)));
            }   
        }
    }  
    return ax2;
}

template <class Type>
matrix<Type> row_operations_result(matrix<Type> const ax, matrix<Type> const b) {
    auto b2 = b;
    
    for (uint64_t k = 0; k < ax.m() - 1; k++) {
        for (uint64_t i = k + 1; i < ax.m(); i++) {
            b2.set(i, 0, b2.get(i, 0) - (ax.get(i, k) * b2.get(k, 0)));
        }
    }
    return b2;
}


template <class Type>
matrix<Type> back_substitution_gauss(matrix<Type> const ax, matrix<Type> const b) {
    auto b2 = b;

    uint64_t m = ax.m();
    uint64_t n = ax.n();
    
    b2.set(n - 1, 0, b2.get(n - 1, 0) / ax.get(n - 1, n - 1));

    for (uint64_t i = m - 2; i < m - 1; i--) {
        Type s = b2.get(i, 0);
        for (uint64_t j = i + 1; j < m; j++) {
            s -= (ax.get(i, j) * b2.get(j, 0));
        }
        b2.set(i, 0, s / ax.get(i, i));    
    }   

    return b2;
}


template <class Type>
matrix<Type> naive_gauss_solve(matrix<Type> const ax, matrix<Type> const b) {
    auto ax2 = naive_gauss_elimination(ax);
    auto b2 = row_operations_result(ax2, b);
    return back_substitution_gauss(ax2, b2);
}


#endif
