#pragma once

#include "cg_base/ijk.hpp"
#include "cg_base/ijkpair.hpp"
#include <cstddef>
#include <string>
#include <vector>

class ijkPairM {
public:
    ijk a, b;
    std::size_t m;

    ijkPairM();
    ijkPairM(ijk a, ijk b, std::size_t m);

    std::string to_string() const;
    std::size_t single_pair_size() const;
    std::size_t single_index() const;

    std::size_t cnt_zero() const;
    std::size_t cnt_nonz() const;

    std::size_t  operator[](std::size_t pos) const;
    std::size_t& operator[](std::size_t pos);

    std::vector<std::size_t> nonzero_index() const;

    std::vector<ijkPairM> hgp_vrr_need(std::size_t pos) const;
};

bool operator< (const ijkPairM &a, const ijkPairM &b);
bool operator> (const ijkPairM &a, const ijkPairM &b);
bool operator<=(const ijkPairM &a, const ijkPairM &b);
bool operator>=(const ijkPairM &a, const ijkPairM &b);
bool operator==(const ijkPairM &a, const ijkPairM &b);
bool operator!=(const ijkPairM &a, const ijkPairM &b);

// hash function of class ijk
struct ijkPairM_hash {
    std::size_t operator()(const ijkPairM &pm) const {
        return hash_val(pm.a.i, pm.a.j, pm.a.k,
                        pm.b.i, pm.b.j, pm.b.k, pm.m);
    }
};

// generate all ijk that i+j+k == n
std::vector<ijkPairM> generate_ijkpair0(std::size_t na, std::size_t nb);