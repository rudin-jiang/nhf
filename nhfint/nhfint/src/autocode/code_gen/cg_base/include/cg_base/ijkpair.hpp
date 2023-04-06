#pragma once

#include "cg_base/hash.hpp"
#include "cg_base/ijk.hpp"
#include <cstddef>
#include <vector>
#include <string>

class ijkPair {
public:
    ijk a, b;
    
    ijkPair() {}
    ijkPair(ijk a, ijk b): a(a), b(b) {}

    std::size_t sum() const { return a.sum() + b.sum(); }

    std::string to_string() const;
    std::string to_comment() const;
    std::size_t pair_size() const;
    std::size_t index() const;

    // form in axbx[idx]
    std::string array_element(std::size_t indexLen) const;

    std::size_t  operator[](std::size_t pos) const;
    std::size_t& operator[](std::size_t pos);

    // hgp hrr
    std::vector<ijkPair> hgp_hrr_need(std::size_t pos) const;
};

bool operator< (const ijkPair &ab, const ijkPair &cd);
bool operator> (const ijkPair &ab, const ijkPair &cd);
bool operator<=(const ijkPair &ab, const ijkPair &cd);
bool operator>=(const ijkPair &ab, const ijkPair &cd);
bool operator==(const ijkPair &ab, const ijkPair &cd);
bool operator!=(const ijkPair &ab, const ijkPair &cd);

// hash function of class ijkPair
struct ijkPair_hash {
    std::size_t operator()(const ijkPair &p) const {
        return hash_val(p.a.i, p.a.j, p.a.k,
                        p.b.i, p.b.j, p.b.k);
    }
};

std::vector<ijkPair> generate_pair(std::size_t na, std::size_t nb);