#include "nhfint/int2e_class_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int2eClassCartD0::Int2eClassCartD0()
: ang0(0), ang1(0), ang2(0), ang3(0),
  nbs0(0), nbs1(0), nbs2(0), nbs3(0),
  num0(0), num1(0), num2(0), num3(0) {}

Int2eClassCartD0::Int2eClassCartD0(
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

    e2D0 = std::vector<double>(num0);
}

double  Int2eClassCartD0::d0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) const {
    assert(i < nbs0); assert(j < nbs1);
    assert(k < nbs2); assert(l < nbs3);
    return e2D0[i * num1 + j * num2 + k * num3 + l];
}

double& Int2eClassCartD0::d0(
    std::size_t i, std::size_t j, 
    std::size_t k, std::size_t l
) {
    assert(i < nbs0); assert(j < nbs1);
    assert(k < nbs2); assert(l < nbs3);
    return e2D0[i * num1 + j * num2 + k * num3 + l];
}

// TODO: two index and one index imp

} // namespace (nhfInt)
