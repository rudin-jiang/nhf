#include "nhfint/int1e_tensor.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int1eTensorD1::Int1eTensorD1(): nBasis(0), pageSize(0) {}

Int1eTensorD1::Int1eTensorD1(std::size_t nBasis)
: nBasis(nBasis), pageSize(nBasis * nBasis) {
    e1D1x = std::vector<double>(pageSize);
    e1D1y = std::vector<double>(pageSize);
    e1D1z = std::vector<double>(pageSize);
}

// derivative of first basis function
double  Int1eTensorD1::d1x0(std::size_t i, std::size_t j) const {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1x[i * nBasis + j];
}

double& Int1eTensorD1::d1x0(std::size_t i, std::size_t j) {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1x[i * nBasis + j];
}

double  Int1eTensorD1::d1y0(std::size_t i, std::size_t j) const {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1y[i * nBasis + j];
}

double& Int1eTensorD1::d1y0(std::size_t i, std::size_t j) {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1y[i * nBasis + j];
}

double  Int1eTensorD1::d1z0(std::size_t i, std::size_t j) const {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1z[i * nBasis + j];
}

double& Int1eTensorD1::d1z0(std::size_t i, std::size_t j) {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1z[i * nBasis + j];
}


// specify the basis function to be derived
double  Int1eTensorD1::d1x(std::size_t a, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1x[(a == 0) ? (i * nBasis + j)
                          : (j * nBasis + i)];
}

double& Int1eTensorD1::d1x(std::size_t a, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1x[(a == 0) ? (i * nBasis + j)
                          : (j * nBasis + i)];
}

double  Int1eTensorD1::d1y(std::size_t a, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1y[(a == 0) ? (i * nBasis + j)
                          : (j * nBasis + i)];
}

double& Int1eTensorD1::d1y(std::size_t a, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1y[(a == 0) ? (i * nBasis + j)
                          : (j * nBasis + i)];
}

double  Int1eTensorD1::d1z(std::size_t a, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1z[(a == 0) ? (i * nBasis + j)
                          : (j * nBasis + i)];
}

double& Int1eTensorD1::d1z(std::size_t a, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D1z[(a == 0) ? (i * nBasis + j)
                          : (j * nBasis + i)];
}


} // namespace (nhfInt)
