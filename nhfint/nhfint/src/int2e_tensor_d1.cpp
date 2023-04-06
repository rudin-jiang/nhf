#include "nhfint/int2e_tensor.hpp"
#include "nhfint/utility.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int2eTensorD1::Int2eTensorD1()
: nBasis(0), pageSize(0), cubeSize(0), pileSize(0) {}

Int2eTensorD1::Int2eTensorD1(std::size_t nBasis): nBasis(nBasis) {
    assert(nBasis > 0);
    pageSize = idx2(nBasis-1, nBasis-1) + 1;
    cubeSize = nBasis * pageSize;
    pileSize = nBasis * cubeSize;

    e2D1x = std::vector<double>(pileSize);
    e2D1y = std::vector<double>(pileSize);
    e2D1z = std::vector<double>(pileSize);
}

// derivative of first basis function
double  Int2eTensorD1::d1x0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) const {
    assert(i < nBasis); assert(j < nBasis); 
    assert(k < nBasis); assert(l < nBasis);
    return e2D1x[i * cubeSize + j * pageSize + idx2(k, l)];
}

double& Int2eTensorD1::d1x0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) {
    assert(i < nBasis); assert(j < nBasis); 
    assert(k < nBasis); assert(l < nBasis);
    return e2D1x[i * cubeSize + j * pageSize + idx2(k, l)];
}

double  Int2eTensorD1::d1y0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) const {
    assert(i < nBasis); assert(j < nBasis); 
    assert(k < nBasis); assert(l < nBasis);
    return e2D1y[i * cubeSize + j * pageSize + idx2(k, l)];
}

double& Int2eTensorD1::d1y0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) {
    assert(i < nBasis); assert(j < nBasis); 
    assert(k < nBasis); assert(l < nBasis);
    return e2D1y[i * cubeSize + j * pageSize + idx2(k, l)];
}

double  Int2eTensorD1::d1z0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) const {
    assert(i < nBasis); assert(j < nBasis); 
    assert(k < nBasis); assert(l < nBasis);
    return e2D1z[i * cubeSize + j * pageSize + idx2(k, l)];
}

double& Int2eTensorD1::d1z0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) {
    assert(i < nBasis); assert(j < nBasis); 
    assert(k < nBasis); assert(l < nBasis);
    return e2D1z[i * cubeSize + j * pageSize + idx2(k, l)];
}



// a = 0:  idx = i * cubeSize + j * pageSize + idx2(k, l)
// a = 1:  idx = j * cubeSize + i * pageSize + idx2(k, l)
// a = 2:  idx = k * cubeSize + l * pageSize + idx2(i, j)
// a = 3:  idx = l * cubeSize + k * pageSize + idx2(i, j)
std::size_t Int2eTensorD1::_d1_index(
    std::size_t a,
    std::size_t i, std::size_t j,
    std::size_t k, std::size_t l
) const {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);

    // TODO: may be using some fast index
    if (a == 0) return i * cubeSize + j * pageSize + idx2(k, l);
    if (a == 1) return j * cubeSize + i * pageSize + idx2(k, l);
    if (a == 2) return k * cubeSize + l * pageSize + idx2(i, j);
    return l * cubeSize + k * pageSize + idx2(i, j);
}

// specify the basis function to be derived
double  Int2eTensorD1::d1x(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);
    return e2D1x[_d1_index(a, i, j, k, l)];   
}

double& Int2eTensorD1::d1x(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);
    return e2D1x[_d1_index(a, i, j, k, l)]; 
}

double  Int2eTensorD1::d1y(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);
    return e2D1y[_d1_index(a, i, j, k, l)]; 
}

double& Int2eTensorD1::d1y(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);
    return e2D1y[_d1_index(a, i, j, k, l)]; 
}

double  Int2eTensorD1::d1z(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);
    return e2D1z[_d1_index(a, i, j, k, l)]; 
}

double& Int2eTensorD1::d1z(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(a < 4);
    assert(i < nBasis); assert(j < nBasis);
    assert(k < nBasis); assert(l < nBasis);
    return e2D1z[_d1_index(a, i, j, k, l)]; 
}


} // namespace (nhfInt)
