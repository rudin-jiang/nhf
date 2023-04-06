#pragma once

#include <cstddef>
#include <vector>


class ijkm {
public:
    std::size_t i, j, k;
    std::size_t m;

    ijkm(): i(0), j(0), k(0), m(0) {}

    ijkm(std::size_t i, std::size_t j, 
         std::size_t k, std::size_t m);
};


