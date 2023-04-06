#pragma once

#include <nhfint/basis.hpp>
#include <nhfmath/matrix.hpp>

namespace nhfInt {

namespace tho {
nhfMath::Matrix int1e_ovlp_d0_cart(const SegmentBasisSet &bsSet);
nhfMath::Matrix int1e_ovlp_d1_cart(const SegmentBasisSet &bsSet);
nhfMath::Matrix int1e_ovlp_d2_cart(const SegmentBasisSet &bsSet);
}  // namespace (tho)


namespace mdx {
nhfMath::Matrix int1e_ovlp_d0_cart(const SegmentBasisSet &bsSet);
nhfMath::Matrix int1e_ovlp_d1_cart(const SegmentBasisSet &bsSet);
nhfMath::Matrix int1e_ovlp_d2_cart(const SegmentBasisSet &bsSet);
}  // namespace (tho)


} // namespace (nhfInt)
