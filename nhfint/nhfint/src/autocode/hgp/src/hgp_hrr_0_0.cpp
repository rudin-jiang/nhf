#include "hgp/hgp_hrr.hpp"
#include <vector>

namespace nhfInt {
namespace hgp {

void hgp_hrr_0_0_(std::vector<double> &a0b0, const std::vector<double> &hrrInp, double x, double y, double z)
{
    for (std::size_t i = 0; i < 1; ++i) {
        a0b0[i] = hrrInp[i];
    }
}

} // namespace (hgp)
} // namespace (nhfInt)
