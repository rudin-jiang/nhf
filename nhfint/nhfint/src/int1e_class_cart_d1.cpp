#include "nhfint/int1e_class_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int1eClassCartD1::Int1eClassCartD1()
: ang0(0), ang1(0), nbs0(0), nbs1(0), num0(0), num1(0) {}

Int1eClassCartD1::Int1eClassCartD1(std::size_t ang0, std::size_t ang1)
: ang0(ang0), ang1(ang1)
{
    assert(ang0 <= AngMomCartMaxL);
    assert(ang1 <= AngMomCartMaxL);

    nbs0 = (ang0 + 1) * (ang0 + 2) / 2;
    nbs1 = (ang1 + 1) * (ang1 + 2) / 2;
    
    num1 = nbs1;
    num0 = nbs0 * num1;
    
    e1D1x = std::vector<double>(2 * num0);
    e1D1y = std::vector<double>(2 * num0);
    e1D1z = std::vector<double>(2 * num0);
}


// Derivatives of first basis function
double  Int1eClassCartD1::d1x0(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1x[i * num1 + j];
}

double& Int1eClassCartD1::d1x0(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1x[i * num1 + j];
}

double  Int1eClassCartD1::d1y0(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1y[i * num1 + j];
}

double& Int1eClassCartD1::d1y0(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1y[i * num1 + j];
}

double  Int1eClassCartD1::d1z0(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1z[i * num1 + j];
}

double& Int1eClassCartD1::d1z0(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1z[i * num1 + j];
}


// Derivatives of second basis function
double  Int1eClassCartD1::d1x1(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1x[num0 + i * num1 + j];
}

double& Int1eClassCartD1::d1x1(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1x[num0 + i * num1 + j];
}

double  Int1eClassCartD1::d1y1(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1y[num0 + i * num1 + j];
}

double& Int1eClassCartD1::d1y1(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1y[num0 + i * num1 + j];
}

double  Int1eClassCartD1::d1z1(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1z[num0 + i * num1 + j];
}

double& Int1eClassCartD1::d1z1(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1z[num0 + i * num1 + j];
}


// Specify the basis function to be derived
double  Int1eClassCartD1::d1x(std::size_t a, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1x[a * num0 + i * num1 + j];
}

double& Int1eClassCartD1::d1x(std::size_t a, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1x[a * num0 + i * num1 + j];
}

double  Int1eClassCartD1::d1y(std::size_t a, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1y[a * num0 + i * num1 + j];
}

double& Int1eClassCartD1::d1y(std::size_t a, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1y[a * num0 + i * num1 + j];
}

double  Int1eClassCartD1::d1z(std::size_t a, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1z[a * num0 + i * num1 + j];
}

double& Int1eClassCartD1::d1z(std::size_t a, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D1z[a * num0 + i * num1 + j];
}


} // namespace (nhfInt)
