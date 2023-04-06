#include "nhfint/int2e_class_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int2eClassCartD2::Int2eClassCartD2()
: ang0(0), ang1(0), ang2(0), ang3(0),
  nbs0(0), nbs1(0), nbs2(0), nbs3(0),
  num0(0), num1(0), num2(0), num3(0) {}

Int2eClassCartD2::Int2eClassCartD2(
    std::size_t ang0, std::size_t ang1,
    std::size_t ang2, std::size_t ang3
) : ang0(0), ang1(0), ang2(0), ang3(0) {
    assert(ang0 <= 2 * AngMomCart);
    assert(ang1 <= 2 * AngMomCart);
    assert(ang2 <= 2 * AngMomCart);
    assert(ang3 <= 2 * AngMomCart);

    nbs0 = (ang0 + 1) * (ang0 + 2) / 2;
    nbs1 = (ang1 + 1) * (ang1 + 2) / 2;
    nbs2 = (ang2 + 1) * (ang2 + 2) / 2;
    nbs3 = (ang3 + 1) * (ang3 + 2) / 2;
    
    num3 = nbs3;
    num2 = nbs2 * num3;
    num1 = nbs1 * num2;
    num0 = nbs0 * num1;

    e2D2xx = std::vector<double>(16 * num0);
    e2D2yy = std::vector<double>(16 * num0);
    e2D2zz = std::vector<double>(16 * num0);
    e2D2xy = std::vector<double>(16 * num0);
    e2D2yz = std::vector<double>(16 * num0);
    e2D2zx = std::vector<double>(16 * num0);
}


// Specify the basis function to be derived
double  Int2eClassCartD2::d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2xx[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2xx[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2yy[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2yy[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2zz[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2zz[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2xy[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2xy[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2yz[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2yz[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2zx[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    return e2D2zx[(4*a + b) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    // d2xy(b, a, i, j, k, l)
    return e2D2xy[(4*b + a) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    // d2xy(b, a, i, j, k, l)
    return e2D2xy[(4*b + a) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    // d2yz(b, a, i, j, k, l)
    return e2D2yz[(4*b + a) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    // d2yz(b, a, i, j, k, l)
    return e2D2yz[(4*b + a) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD2::d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    // d2zx(b, a, i, j, k, l)
    return e2D2zx[(4*b + a) * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD2::d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    assert(a < 4); assert(b < 4); assert(i < nbs0); assert(j < nbs1); assert(k < nbs2); assert(l < nbs3);
    // d2zx(b, a, i, j, k, l)
    return e2D2zx[(4*b + a) * num0 + i * num1 + j * num2 + k * num3 + l];
}


} // namespace (nhfInt)
