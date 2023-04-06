#include "nhfmath/matrix.hpp"
#include <cstddef>
#include <nhfint/int1e_ovlp_cart.hpp>


namespace nhfInt {

namespace tho {

nhfMath::Matrix int1e_ovlp_cart(const SegmentBasisSet &bsSet) {
    const std::size_t nBasis = bsSet.cart_basis_size();
    const std::size_t nShell = bsSet.shell_size();
    nhfMath::Matrix ovlp(nBasis, nBasis);


    // TODO:

    // for (std::size_t i = 0; i < nShell; ++i) {
    // for (std::size_t j = i; j < nShell; ++j) {
        
    //     for (std::size_t ia = 0; ia < bsSet[i].basis_size(); ++ia) {
    //     for (std::size_t jb = 0; jb < bsSet[j].basis_size(); ++jb) {


    //     }}

    // }}

    return ovlp;
}

}  // namespace (nhfInt)


} // namespace (nhfInt)