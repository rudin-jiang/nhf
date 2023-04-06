#include "nhfint/int2e_class_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int2eClassCartD1::Int2eClassCartD1()
: ang0(0), ang1(0), ang2(0), ang3(0),
  nbs0(0), nbs1(0), nbs2(0), nbs3(0),
  num0(0), num1(0), num2(0), num3(0) {}

Int2eClassCartD1::Int2eClassCartD1(
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

    e2D1x = std::vector<double>(4 * num0);
    e2D1y = std::vector<double>(4 * num0);
    e2D1z = std::vector<double>(4 * num0);
}


double  Int2eClassCartD1::d1x0(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1x0(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1y0(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1y0(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1z0(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1z0(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1x1(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1x1(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1y1(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1y1(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1z1(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1z1(
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1x2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[2 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1x2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[2 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1y2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[2 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1y2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[2 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1z2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[2 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1z2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[2 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1x3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[3 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1x3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[3 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1y3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[3 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1y3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[3 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1z3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[3 * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1z3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[3 * num0 + i * num1 + j * num2 + k * num3 + l];
}


// Specify the basis function to be derived
double  Int2eClassCartD1::d1x(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(a < 4); 
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[a * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1x(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) {
    assert(a < 4); 
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1x[a * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1y(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const {
    assert(a < 4); 
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[a * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1y(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(a < 4); 
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1y[a * num0 + i * num1 + j * num2 + k * num3 + l];
}

double  Int2eClassCartD1::d1z(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) const {
    assert(a < 4); 
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[a * num0 + i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD1::d1z(
    std::size_t a, 
    std::size_t i, std::size_t j, 
    std::size_t k ,std::size_t l
) {
    assert(a < 4); 
    assert(i < nbs0); assert(j < nbs1); 
    assert(k < nbs2); assert(l < nbs3);
    return e2D1z[a * num0 + i * num1 + j * num2 + k * num3 + l];
}


} // namespace (nhfInt)
