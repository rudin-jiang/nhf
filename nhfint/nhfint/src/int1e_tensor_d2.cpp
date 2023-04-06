#include "nhfint/int1e_tensor.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int1eTensorD2::Int1eTensorD2(): nBasis(0), pageSize(0) {}

Int1eTensorD2::Int1eTensorD2(std::size_t nBasis)
: nBasis(nBasis), pageSize(nBasis * nBasis) {
    e1D2xx = std::vector<double>(2 * pageSize);
    e1D2yy = std::vector<double>(2 * pageSize);
    e1D2zz = std::vector<double>(2 * pageSize);
    e1D2xy = std::vector<double>(2 * pageSize);
    e1D2yz = std::vector<double>(2 * pageSize);
    e1D2zx = std::vector<double>(2 * pageSize);
}

// a = 0, b = 0 : [i * nBasis + j]
// a = 0, b = 1 : [pageSize + i * nBasis + j]
// a = 1, b = 0 : [pageSize + j * nBasis + i]
// a = 1, b = 1 : [j * pageSize + i]
std::size_t Int1eTensorD2::_d2_index(
    std::size_t a, std::size_t b,
    std::size_t i, std::size_t j
) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    if (a == 0 && b == 0)   return i * nBasis + j;
    if (a == 0 && b == 1)   return pageSize + i * nBasis + j;
    if (a == 1 && b == 0)   return pageSize + j * nBasis + i;
    /*  a == 1 && b == 1 */ return j * nBasis + i;
}

double  Int1eTensorD2::d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2xx[_d2_index(a, b, i, j)];
}

double& Int1eTensorD2::d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2xx[_d2_index(a, b, i, j)];
}

double  Int1eTensorD2::d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2yy[_d2_index(a, b, i, j)];
}

double& Int1eTensorD2::d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2yy[_d2_index(a, b, i, j)];
}

double  Int1eTensorD2::d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2zz[_d2_index(a, b, i, j)];
}

double& Int1eTensorD2::d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2zz[_d2_index(a, b, i, j)];
}

double  Int1eTensorD2::d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2xy[_d2_index(a, b, i, j)];
}

double& Int1eTensorD2::d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2xy[_d2_index(a, b, i, j)];
}

double  Int1eTensorD2::d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2yz[_d2_index(a, b, i, j)];
}

double& Int1eTensorD2::d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2yz[_d2_index(a, b, i, j)];
}

double  Int1eTensorD2::d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2zx[_d2_index(a, b, i, j)];
}

double& Int1eTensorD2::d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    return e1D2zx[_d2_index(a, b, i, j)];
}

double  Int1eTensorD2::d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    // d2xy(b, a, i, j)
    return e1D2xy[_d2_index(b, a, i, j)];
}

double& Int1eTensorD2::d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    // d2xy(b, a, i, j)
    return e1D2xy[_d2_index(b, a, i, j)];
}

double  Int1eTensorD2::d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    // d2yz(b, a, i, j)
    return e1D2yz[_d2_index(b, a, i, j)];
}

double& Int1eTensorD2::d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    // d2yz(b, a, i, j)
    return e1D2yz[_d2_index(b, a, i, j)];
}

double  Int1eTensorD2::d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    // d2zx(b, a, i, j)
    return e1D2zx[_d2_index(b, a, i, j)];
}

double& Int1eTensorD2::d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2); assert(b < 2);
    assert(i < nBasis); assert(j < nBasis);
    // d2zx(b, a, i, j)
    return e1D2zx[_d2_index(b, a, i, j)];
}



} // namespace
