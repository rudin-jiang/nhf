#pragma once

#include <cstddef>
#include <vector>
#include <string>

namespace nhfInt {

// When we use a one-dimensional array to store a symmetric matrix, 
// idx2 calculates the position of the matrix element in the array.
inline std::size_t idx2(std::size_t i, std::size_t j)
{ return  i>j ? i * (i+1) / 2 + j : j * (j+1) / 2 + i; }

inline std::size_t idx4(std::size_t i, std::size_t j, 
                        std::size_t k, std::size_t l)
{ return  idx2(idx2(i,j), idx2(k,l)); }


// Normalization constant for Cartesian Gaussian function.
// ref: doi: 10.1143/JPSJ.21.2313 (eq 2.2)
// Overloads for 1D and 3D.
double gauss_norm_cart(double alpha, std::size_t l);

double gauss_norm_cart(double alpha, std::size_t l, 
                                     std::size_t m, 
                                     std::size_t n);


}   // namespace (nhfInt)