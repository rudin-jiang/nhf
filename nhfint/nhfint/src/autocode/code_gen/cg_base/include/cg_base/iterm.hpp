#pragma once

#include "cg_base/hash.hpp"
#include <cstddef>
#include <string>
#include <cassert>
#include <vector>


class ITerm {
public:
    std::size_t i, j, k, l;

    ITerm();
    ITerm(std::size_t i, std::size_t j,
          std::size_t k, std::size_t l);

    std::string to_string() const;
};

bool operator< (const ITerm &a, const ITerm &b);
bool operator> (const ITerm &a, const ITerm &b);
bool operator<=(const ITerm &a, const ITerm &b);
bool operator>=(const ITerm &a, const ITerm &b);
bool operator==(const ITerm &a, const ITerm &b);
bool operator!=(const ITerm &a, const ITerm &b);

struct iterm_hash {
    std::size_t operator()(const ITerm &a) const {
        return hash_val(a.i, a.j, a.k, a.l);
    }
};


// generate all ITerm that
// i range from 0 to na,  j range from 0 to nb,
// k range from 0 to nc,  l range from 0 to nd.
std::vector<ITerm> generate_iterm(std::size_t na, std::size_t nb, 
                                  std::size_t nc, std::size_t nd);