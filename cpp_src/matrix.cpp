#include "matrix.h"
#include "matrix_factorisations/gram_schmidt.h"
#include "matrix_factorisations/modified_gram_schmidt.h"
#include "matrix_factorisations/pivoting_gauss_elimination.h"

#include <complex>

int main(int agrc, char *argv[]) {
    double myints[] = {1,2,3,4,5,6,7,8,9,10,11,12};
    double myints2[] = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    std::vector<double> v(myints, myints + sizeof(myints) / sizeof(double));  
    std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));  

    auto m = matrix<double>(v, 3);
    auto n = matrix<double>(v2, 3);
   
    matrix<double> o = m;
    o.rbind(m);

    matrix<double> p = m;
    p.cbind(m);

    std::cout << (m == n) << std::endl;
    std::cout << (m != n) << std::endl;
     
    m.print();
    std::cout << std::endl;
    m.transpose();
    m.print();
    std::cout << std::endl;

    std::cout << (m == n) << std::endl;
    std::cout << (m != n) << std::endl;

    matrix<double> a = m + m;
    a.print();
    std::cout << std::endl;

    matrix<double> b = m - m;
    b.print();
    std::cout << std::endl;

    matrix<double> c = n * m;
    c.print();
    std::cout << std::endl;

    o.print();
    std::cout << std::endl;
  
    p.print();
    std::cout << std::endl;

    auto q = 5.0 * m;
    q.print();
    std::cout << std::endl;

    m.print();
    std::cout << std::endl;

    
    double myints3[] = {1,1,0,1,0,1,0,1,1};    
    std::vector<double> v3(myints3, myints3 + sizeof(myints3) / sizeof(double));  
    auto z = matrix<double>(v3, 3);
   
    z.print();
    std::cout << std::endl;
    
    auto gs = gram_schmidt(z);
    std::get<0>(gs).print();
    std::cout << std::endl;
    std::get<1>(gs).print();
    std::cout << std::endl;

    
    auto gs2 = modified_gram_schmidt(z);
    std::get<0>(gs2).print();
    std::cout << std::endl;
    std::get<1>(gs2).print();
    std::cout << std::endl;

    return 0;
}
