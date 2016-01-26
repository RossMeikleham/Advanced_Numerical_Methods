#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "../matrix.h"
#include <iostream>
#include <vector>
#include <tuple>


template <class Type>
matrix<Type> cholesky(matrix<Type> a) {
    matrix<Type> l = zeros(a.m(), a.n());
     
    for (uint64_t i = 0; i < a.n(); i++) {
        for (uint64_t j = i; j < a.m(); j++) {
            if (i == j) {
                auto s = a.get(i, i);
                for (uint64_t k = 0; i && k < i; k++) {
                    s -= std::pow(l.get(i, k), 2);  
                }
                l.set(i, i, std::sqrt(s));
            }

            else {
                auto s = a.get(j, i);
                for (uint64_t k = 0; i && k < i; k++) {
                    s -= (l.get(i, k) * l.get(j, k));
                }
                l.set(j, i, s / l.get(i, i));
            }
        } 
    }

    return l;
}
         
#endif
