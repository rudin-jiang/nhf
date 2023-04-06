#include "nhfint/int1e_class_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {


Int1eClassCartD0::Int1eClassCartD0()
: ang0(0), ang1(0), nbs0(0), nbs1(0), num0(0), num1(0) {}

Int1eClassCartD0::Int1eClassCartD0(std::size_t ang0, std::size_t ang1)
: ang0(ang0), ang1(ang1)
{
    assert(ang0 <= AngMomCartMaxL);
    assert(ang1 <= AngMomCartMaxL);

    nbs0 = (ang0 + 1) * (ang0 + 2) / 2;
    nbs1 = (ang1 + 1) * (ang1 + 2) / 2;
    
    num1 = nbs1;
    num0 = nbs0 * num1;
    
    e1D0 = std::vector<double>(num0);
}

double Int1eClassCartD0::d0(std::size_t i, std::size_t j) const {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D0[i * num1 + j];
}

double& Int1eClassCartD0::d0(std::size_t i, std::size_t j) {
    assert(i < nbs0);
    assert(j < nbs1);
    return e1D0[i * num1 + j];
}


} // namespace (nhfInt)