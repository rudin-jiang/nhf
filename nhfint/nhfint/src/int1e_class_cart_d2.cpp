#include "nhfint/int1e_class_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <cassert>

namespace nhfInt {

Int1eClassCartD2::Int1eClassCartD2()
: ang0(0), ang1(0), nbs0(0), nbs1(0), num0(0), num1(0) {}

Int1eClassCartD2::Int1eClassCartD2(std::size_t ang0, std::size_t ang1)
: ang0(ang0), ang1(ang1)
{
    assert(ang0 <= AngMomCartMaxL);
    assert(ang1 <= AngMomCartMaxL);

    nbs0 = (ang0 + 1) * (ang0 + 2) / 2;
    nbs1 = (ang1 + 1) * (ang1 + 2) / 2;
    
    num1 = nbs1;
    num0 = nbs0 * num1;
    
    e1D2xx = std::vector<double> (4 * num0);
    e1D2yy = std::vector<double> (4 * num0);
    e1D2zz = std::vector<double> (4 * num0);
    e1D2xy = std::vector<double> (4 * num0);
    e1D2yz = std::vector<double> (4 * num0);
    e1D2zx = std::vector<double> (4 * num0);
}


// Specify the basis function to be derived
double  Int1eClassCartD2::d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2xx[(2*a + b) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2xx[(2*a + b) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2yy[(2*a + b) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2yy[(2*a + b) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2zz[(2*a + b) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2zz[(2*a + b) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2xy[(2*a + b) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2xy[(2*a + b) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2yz[(2*a + b) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2yz[(2*a + b) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2zx[(2*a + b) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D2zx[(2*a + b) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);

    // d2xy(b, a, i, j)
    return e1D2xy[(2*b + a) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    
    // d2xy(b, a, i, j)
    return e1D2xy[(2*b + a) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    
    // d2yz(b, a, i, j)
    return e1D2yz[(2*b + a) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);

    // d2yz(b, a, i, j)
    return e1D2yz[(2*b + a) * num0 + i * num1 + j];
}

double  Int1eClassCartD2::d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);

    // d2zx(b, a, i, j)
    return e1D2zx[(2*b + a) * num0 + i * num1 + j];
}

double& Int1eClassCartD2::d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) {
    assert(a < 2);
    assert(b < 2);
    assert(i < nbs0);
    assert(j < nbs1);
    
    // d2zx(b, a, i, j)
    return e1D2zx[(2*b + a) * num0 + i * num1 + j];
}


} // namespace (nhfInt)


