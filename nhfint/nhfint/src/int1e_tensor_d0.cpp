#include "nhfint/int1e_tensor.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {

Int1eTensorD0::Int1eTensorD0(): nBasis(0), pageSize(0) {}

Int1eTensorD0::Int1eTensorD0(std::size_t nBasis)
: nBasis(nBasis), pageSize(nBasis * nBasis) {
    e1D0 = std::vector<double>(pageSize);
}

double  Int1eTensorD0::d0(std::size_t i, std::size_t j) const {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D0[i * nBasis + j];
}

double& Int1eTensorD0::d0(std::size_t i, std::size_t j) {
    assert(i < nBasis);
    assert(j < nBasis);
    return e1D0[i * nBasis + j];
}

} // namespace (nhfInt)
