#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <math.h>
#include <complex>

template <class Type> class matrix {
  private:
    uint64_t m_m, m_n;
    std::vector<Type> elements;

  public:
    matrix(std::vector<Type> v, uint64_t n) {
        elements = v;
        m_n = n;
        m_m = elements.size() / n;
    }
     
    uint64_t m() const {return m_m;} 
    uint64_t n() const {return m_n;} 

    Type get(uint64_t y, uint64_t x) const {
        return elements[(y * m_n) + x];
    }

    void set(uint64_t y, uint64_t x, Type t) {
        elements[(y * m_n) + x] = t;
    }

    matrix<Type> col(uint64_t col) const {
        std::vector<Type> v = std::vector<Type>(m_m, 0);
       
        #pragma omp parallel for 
        for (uint64_t y = 0; y < m_m; y++) {
            v[y] = elements[(y * m_n) + col];          
        }

        return matrix<Type>(v, 1);
    }

    matrix<Type> row(uint64_t row) const {
        std::vector<Type> v = std::vector<Type>(m_n, 0);
        
        #pragma omp parallel for
        for (uint64_t x = 0; x < m_n; x++) {
            v[x] = elements[(row * m_n) + x];
        }

        return matrix<Type>(v, m_n);
    }

    void setRow(uint64_t row, matrix<Type> m) {
        #pragma omp parallel for
        for (uint64_t x = 0; x < m_n; x++) {
            elements[(row * m_n) + x] = m.get(0, x);            
        }
    }

    void transpose() {
        std::vector<Type> temp_v = elements;

        #pragma omp parallel for
        for (uint64_t y = 0; y < m_n; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m_m; x++) {
                elements[(y * m_m) + x] = temp_v[(x * m_n) + y]; 
            }
        } 
          
        uint64_t temp = m_m;
        m_m = m_n;
        m_n = temp;
    }

    bool operator==(matrix<Type> const &m2) const {
        
        // Check dimensional equality
        if (m_m != m2.m_m || m_n != m2.m_n) {
            return false;
        }

        // Check each element is equal
        for (uint64_t y = 0; y < m_m; y++) {
            for (uint64_t x = 0; x < m_n; x++) {
                if (elements[(y * m_n) + x] != m2.get(y, x)) {
                    return false;
                }
            } 
        }

        return true; 
    }

    bool equals(matrix<Type> const &m2, Type delta) {
        // Check dimensional equality
        if (m_m != m2.m_m || m_n != m2.m_n) {
            return false;
        }

        // Check each element is equal
        for (uint64_t y = 0; y < m_m; y++) {
            for (uint64_t x = 0; x < m_n; x++) {
                if (std::abs(elements[(y * m_n) + x] - m2.get(y, x)) > delta) {
                    return false;
                }
            } 
        }
        return true;
    }

    bool operator!=(matrix<Type> const &m2) const {
        return !(*this == m2);
    }

    matrix<Type> operator+(matrix<Type> const &m2) const {
        
        std::vector<Type> new_v = elements;
        #pragma omp parallel for
        for (uint64_t y = 0; y < m_m; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m_n; x++) {
                new_v[(y * m_n) + x] = elements[(y * m_n) + x] + m2.get(y, x);
            }
        }
        
        return matrix(new_v, m_n);    
    }

    void add(matrix<Type> const &m2) {
        #pragma omp parallel for
        for (uint64_t y = 0; y < m_m; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m_n; x++) {
                elements[(y * m_n) + x] += m2.get(y, x);
            }
        }
    }

    matrix<Type> operator-(matrix<Type> const &m2) const {
        
        std::vector<Type> new_v = elements;
        #pragma omp parallel for
        for (uint64_t y = 0; y < m_m; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m_n; x++) {
                new_v[(y * m_n) + x] = elements[(y * m_n) + x] - m2.get(y, x);
            }
        }
        
        return matrix(new_v, m_n);    
    }

    void sub(matrix<Type> const &m2) {
            #pragma omp parallel for
            for (uint64_t y = 0; y < m_m; y++) {
                #pragma omp parallel for
                for (uint64_t x = 0; x < m_n; x++) {
                    elements[(y * m_n) + x] -= m2.get(y, x);
                }
            }
    }

    matrix<Type> operator*(matrix<Type> const &m2) const {
        std::vector<Type> new_v = std::vector<Type>();
        new_v.resize(m_m * m2.m_n);

        #pragma omp parallel for
        for (uint64_t y = 0; y < m_m; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m2.m_n; x++) {
                Type res = 0;
                for (uint64_t o = 0; o < m_n; o++) {
                    res += this->get(y, o) * m2.get(o, x);
                }
                new_v[(y * m2.m_n) + x] = res;
            }
        }
        
        return matrix(new_v, m2.n());    
    }
    
    void mult(matrix<Type> const &m2) {
        std::vector<Type> new_v = std::vector<Type>();
        new_v.resize(m_m * m2.m_n);

        #pragma omp parallel for
        for (uint64_t y = 0; y < m_m; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m2.m_n; x++) {
                Type res = 0;
                for (uint64_t o = 0; o < m_n; o++) {
                    res += this->get(y, o) * m2.get(o, x);
                }
                new_v[(y * m2.n()) + x] = res;
            }
        }
     
        elements = new_v;
        m_n = m2.m_n;   
    }
     
    void scalarMult(Type const &c) {
        #pragma omp parallel for
        for (uint64_t y = 0; y < m_m; y++) {
            #pragma omp parallel for
            for (uint64_t x = 0; x < m_n; x++) {
                elements[(y * m_n) + x] *= c;
            }
        }
    }

    // Given a matrix with the same number of columns,
    // appends this matrix "below" the current matrix to
    // create a new matrix 
    void rbind(matrix<Type> const &m2) {
        uint64_t orig_sz = elements.size();
        uint64_t append_sz = m2.elements.size();
        elements.resize(orig_sz + append_sz);
        
        #pragma omp parallel for
        for (uint64_t i = 0; i < append_sz; i++) {
            elements[orig_sz + i] = m2.elements[i];
        }

        m_m += m2.m_m;
    } 
    
    // Given a matrix with the same number of rows
    // appends this matrix to the right of the current
    // matrix to create a new matrix
    void cbind(matrix<Type> const &m2) {
        matrix<Type> m2_t = m2;
        
        transpose();
        m2_t.transpose();
        
        rbind(m2_t);
        transpose();
    }
    
    // Convert 1x1 matrix to a scalar value
    Type singleton() {
        return elements[0];
    }

    Type abs() {
        Type sum = 0;
        for (uint64_t y = 0; y < m_m; y++) {
            for (uint64_t x = 0; x < m_n; x++) {
                Type val = elements[(y * m_n) + x];
                sum += (val * val);
            }
        }
        return sqrt(sum);
    }

    void print() {
        for (uint64_t y = 0; y < m_m; y++) {
            for (uint64_t x = 0; x < m_n; x++) {
                std::cout << elements[(y * m_n) + x] << " ";
            }
            std::cout << std::endl; 
        }   
    }
};


template <class Type> 
matrix<Type> operator*(Type const &c, matrix<Type> const &m) {
    auto m2 = m;
    m2.scalarMult(c);
    return m2;
}


template <class Type> 
matrix<Type> operator/(matrix<Type> const &m, Type const &c) {
    auto m2 = m;
    for (uint64_t y = 0; y < m.m(); y++) {
        for (uint64_t x = 0; x < m.n(); x++) {
            m2.set(y, x, m2.get(y, x) / c);           
        }
    }
    return m2;
}


matrix<double> zeros(uint64_t m, uint64_t n) {
    std::vector<double> v = std::vector<double>(m * n, 0);   
    return matrix<double>(v, n); 
}


#endif



