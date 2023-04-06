#include "ryspoly/roots_weights.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

namespace nhfInt {
namespace rysPoly {

void rys_roots_weights(
    std::vector<double> &rts,
    std::vector<double> &wts,
    std::size_t n, double x
) {
    assert(rts.size() >= n);
    assert(wts.size() >= n);

    for (std::size_t i = 0; i < n; ++i) {
        rts[i] = wts[i] = 0;
    }
}

} // namespace (rysPoly)
} // namespace (nhfInt)