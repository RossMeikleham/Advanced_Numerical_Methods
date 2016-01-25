#include "../matrix.h"
#include <iostream>
#include <vector>
#include <tuple>


/* gram_schmidt process of factorising a given matrix
 * returns a tuple with the first element containing the "q" matrix
 * and the second element containing the "r" matrix */
template <class Type> 
std::tuple<matrix<Type>, matrix<Type>> modified_gram_schmidt(matrix<Type> a) {

    matrix<Type> r = matrix<Type>(std::vector<Type>(a.n() * a.n(), 0), a.n());
    auto q = std::vector<matrix<Type>>(a.m(), matrix<Type>(std::vector<Type>(1, 0), 1));
    
    for (uint64_t n = 0; n < a.n(); n++) {
       auto an = a.col(n);
       auto v = an;
        
       r.set(n, n, v.abs());
       q[n] = v / r.get(n, n);
        
       // Calc R-Values for n-th column 
       for (uint64_t i = 0; i < n; i++) {
           auto temp = q[i];
           temp.transpose();
           r.set(i, n, (temp * v).singleton());
           v.sub(r.get(i, n) * q[i]);
       } 
    }

    for (uint64_t i = 1; i < a.n(); i++) {
            q[0].cbind(q[i]);
    } 

    return std::make_tuple(q[0], r);
}
