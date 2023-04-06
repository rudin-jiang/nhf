#include "nhfint/int2e_tensor.hpp"
#include "nhfint/utility.hpp"
#include <cstddef>
#include <vector>
#include <cassert>


namespace nhfInt {


Int2eTensorD0::Int2eTensorD0(): nBasis(0), pageSize(0) {}

Int2eTensorD0::Int2eTensorD0(std::size_t nBasis)
: nBasis(nBasis), pageSize(idx2(nBasis-1, nBasis-1) + 1) {
    assert(nBasis > 0);
    std::size_t nInt2eD0 = idx2(pageSize-1, pageSize-1) + 1;
    e2D0 = std::vector<double>(nInt2eD0);
}

double  Int2eTensorD0::d0(std::size_t i, std::size_t j, 
                          std::size_t k, std::size_t l) const {
    assert(i < nBasis);
    assert(j < nBasis);
    assert(k < nBasis);
    assert(l < nBasis);
    return e2D0[idx4(i, j, k, l)];
}

double& Int2eTensorD0::d0(std::size_t i, std::size_t j, 
                          std::size_t k, std::size_t l) {
    assert(i < nBasis);
    assert(j < nBasis);
    assert(k < nBasis);
    assert(l < nBasis);
    return e2D0[idx4(i, j, k, l)];
}


} // namespace (nhfInt)

