#pragma once

#include <cstddef>
#include <vector>

namespace nhfInt {
namespace rysPoly {

// TODO: Implement
void rys_roots_weights(
    std::vector<double> &rts,
    std::vector<double> &wts,
    std::size_t n, double x
);

} // namespace (rysPoly)
} // namespace (nhfInt)