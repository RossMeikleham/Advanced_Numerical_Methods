#define CATCH_CONFIG_MAIN

#include "../include/catch.hpp"
#include "../cpp_src/matrix.h"
#include "../cpp_src/matrix_factorisations/gram_schmidt.h"
#include "../cpp_src/matrix_factorisations/gaussian_elimination.h"
#include "../cpp_src/matrix_factorisations/pivoting_gauss_elimination.h"

#include <vector>

TEST_CASE("Extracting Columns", "[matrix]") {
    double myints[] = {1,2,3,4,5,6,7,8,9};
    std::vector<double> v(myints, myints + sizeof(myints) / sizeof(double));  
    auto m = matrix<double>(v, 3);
   
    SECTION("Extracting middle column from 3x3 matrix") {
        double myints2[] = {2, 5, 8};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 1);

        REQUIRE(m.col(1) == n);    
    }

    SECTION("Extracting Most Left Column") {
        double myints2[] = {1, 4, 7};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 1);

        REQUIRE(m.col(0) == n);    
    } 
   
    SECTION("Extracting Right Most Column") {
        double myints2[] = {3, 6, 9};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 1);

        REQUIRE(m.col(2) == n);    
    } 
}


TEST_CASE("Extracting Rows", "[matrix]") {
    double myints[] = {1,2,3,4,5,6,7,8,9};
    std::vector<double> v(myints, myints + sizeof(myints) / sizeof(double));  
    auto m = matrix<double>(v, 3);
   
    SECTION("Extracting middle row from 3x3 matrix") {
        double myints2[] = {4, 5, 6};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        REQUIRE(m.row(1) == n);    
    }

    SECTION("Extracting Top Row") {
        double myints2[] = {1, 2, 3};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        REQUIRE(m.row(0) == n);    
    } 
   
    SECTION("Extracting Bottom Row") {
        double myints2[] = {7, 8, 9};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        REQUIRE(m.row(2) == n);    
    } 
}


TEST_CASE("Adding", "[matrix]") {
    double myints[] = {1,2,3,4,5,6,7,8,9};
    std::vector<double> v(myints, myints + sizeof(myints) / sizeof(double));  
    auto m = matrix<double>(v, 3);
   
    SECTION("Identity Adding") {
        double myints2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        REQUIRE(m == m + n);    
    }
    
    SECTION("Identity Adding 2") {
        double myints2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        auto m2 = m;
        m.add(n);
        REQUIRE(m == m2);    
    }
    
    SECTION("Adding Numbers") {
        double myints2[] = {23, 44, 92, 4324, 6536, 4312, 212, 213431, 87588};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        for (uint64_t i = 0; i < 9; i++) {
            v[i] += v2[i];  
        }

        auto sum_m = matrix<double>(v, 3);
        REQUIRE(m + n == sum_m); 
    
    }

    SECTION("Adding Numbers 2") {
        double myints2[] = {23, 44, 92, 4324, 6536, 4312, 212, 213431, 87588};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        for (uint64_t i = 0; i < 9; i++) {
            v[i] += v2[i];  
        }

        auto sum_m = matrix<double>(v, 3);
        m.add(n);
        REQUIRE(m == sum_m); 
    }

}


TEST_CASE("Subtracting", "[matrix]") {
    double myints[] = {1,2,3,4,5,6,7,8,9};
    std::vector<double> v(myints, myints + sizeof(myints) / sizeof(double));  
    auto m = matrix<double>(v, 3);
   
    SECTION("Identity Subtracting") {
        double myints2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        REQUIRE(m == m - n);    
    }
    
    SECTION("Identity Subtracting 2") {
        double myints2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        auto m2 = m;
        m.sub(n);
        REQUIRE(m == m2);    
    }
    
    SECTION("Subtracting Numbers") {
        double myints2[] = {39749832, 921374, 4, 3276, 2789321, 4289, 2304, 342, 742};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        for (uint64_t i = 0; i < 9; i++) {
            v[i] -= v2[i];  
        }

        auto sum_m = matrix<double>(v, 3);
        REQUIRE(m - n == sum_m); 
    
    }

    SECTION("Adding Numbers 2") {
        double myints2[] = {39749832, 921374, 4, 3276, 2789321, 4289, 2304, 342, 742};   
        std::vector<double> v2(myints2, myints2 + sizeof(myints2) / sizeof(double));
        auto n = matrix<double>(v2, 3);

        for (uint64_t i = 0; i < 9; i++) {
            v[i] -= v2[i];  
        }

        auto sum_m = matrix<double>(v, 3);
        m.sub(n);
        REQUIRE(m == sum_m); 
    }

}



TEST_CASE("Naive Gauss Elimination", "[matrix]") {
    double ax_a[] = {1, 5, -2, -7};
    double b_a[] = {7, -5};

    std::vector<double> ax_v(ax_a, ax_a + sizeof(ax_a) / sizeof(double));  
    std::vector<double> b_v(b_a, b_a + sizeof(b_a) / sizeof(double));  
    
    auto ax = matrix<double>(ax_v, 2);
    auto b = matrix<double>(b_v, 1);
  
    SECTION("Naive Gauss Elim") {
        double soln_a[] = {1, 5, -2, 3};
        auto soln = matrix<double>(std::vector<double>(soln_a, soln_a + sizeof(soln_a) / sizeof(double)), 2);
        auto result = naive_gauss_elimination(ax);
       
        REQUIRE(result == soln);

    }
 
    SECTION("Naive Gauss Solving") {
        double soln_a[] = {-8, 3};
        auto soln = matrix<double>(std::vector<double>(soln_a, soln_a + sizeof(soln_a) / sizeof(double)), 1); 
        auto result = naive_gauss_solve(ax, b);
        
        REQUIRE(result == soln);    
    }
}

/*  
TEST_CASE("Partial Pivoting Gauss Elimination", "[matrix]") {
    SECTION("2x3") {    
        double ax_a[] = {4, 2, 5, 5, 6, 1};
        std::vector<double> ax_v(ax_a, ax_a + sizeof(ax_a) / sizeof(double));  
        auto ax = matrix<double>(ax_v, 3);
    
        double actual_a[] = {5, 6, 1, 0, -2.8, 4.2};
        std::vector<double> actual_v(actual_a, actual_a + sizeof(actual_a) / sizeof(double));
        auto actual = matrix<double>(actual_v, 3);
         
        REQUIRE(pivoting_gauss_elimination(ax) == actual);    
    }
*/

/*  
    SECTION("3x4") {
        double ax_a[] = {3, 2, -4, 3, 2, 3, 3, 15, 5, -3, 1, 14};
    }
    SECTION("4x5") { 
        double ax_a[] = {2, 1, 1, 0, 1, 4, 3, 3, 1, 2, 8, 7, 9, 5, 3, 6, 7, 9, 8, 4};
    }
*/
