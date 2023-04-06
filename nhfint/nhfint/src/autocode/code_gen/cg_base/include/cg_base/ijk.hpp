#pragma once

#include "cg_base/hash.hpp"
#include <cassert>
#include <cstddef>
#include <regex>
#include <vector>
#include <string>

class ijk {
public:
    std::size_t i, j, k;

    ijk(): i(0), j(0), k(0) {}
    ijk(std::size_t i, std::size_t j, std::size_t k);
    ijk(std::string line);

    std::size_t sum() const;
    std::size_t cnt_zero() const;
    std::size_t cnt_nonz() const;

    std::vector<std::size_t> nonzero_index() const;
    std::size_t first_nonzero_index() const;

    std::string to_string() const;
    std::string to_comment() const;
    std::size_t shell_size() const;
    std::size_t index() const;

    std::size_t  operator[](std::size_t pos) const;
    std::size_t& operator[](std::size_t pos);
};

bool operator< (const ijk &a, const ijk &b);
bool operator> (const ijk &a, const ijk &b);
bool operator<=(const ijk &a, const ijk &b);
bool operator>=(const ijk &a, const ijk &b);
bool operator==(const ijk &a, const ijk &b);
bool operator!=(const ijk &a, const ijk &b);

// hash function of class ijk
struct ijk_hash {
    std::size_t operator()(const ijk &a) const {
        return hash_val(a.i, a.j, a.k);
    }
};

// generate all ijk that i+j+k == n
std::vector<ijk> generate_shell(std::size_t n);



