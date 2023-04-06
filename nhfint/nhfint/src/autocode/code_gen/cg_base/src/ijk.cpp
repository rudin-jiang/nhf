#include "cg_base/ijk.hpp"
#include "cg_base/utility.hpp"
#include <cstddef>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <vector>
#include <set>


ijk::ijk(std::size_t i, std::size_t j, std::size_t k)
: i(i), j(j), k(k) {
    assert(i < 100);
    assert(j < 100);
    assert(k < 100);
}

ijk::ijk(std::string line) {
    std::istringstream iss(line);
    iss >> i >> j >> k;

    assert(i < 100);
    assert(j < 100);
    assert(k < 100);
}

std::size_t ijk::sum() const { return i + j + k; }
std::size_t ijk::cnt_zero() const { return !i + !j + !k; }
std::size_t ijk::cnt_nonz() const { return !!i + !!j + !!k; }

std::vector<std::size_t> ijk::nonzero_index() const {
    std::vector<std::size_t> ret;
    if (i != 0) ret.push_back(0);
    if (j != 0) ret.push_back(1);
    if (k != 0) ret.push_back(2);
    return ret;
}

std::size_t ijk::first_nonzero_index() const {
    assert(cnt_nonz() > 0);
    if (i != 0) return 0;
    if (j != 0) return 1;
    if (k != 0) return 2;
    return -1;
}

std::string ijk::to_string() const {
    char buff[64];
    std::sprintf(buff, "%2zu  %2zu  %2zu", i, j, k);
    return std::string(buff);
}

std::string ijk::to_comment() const {
    return "(" + encode_number(i) 
               + encode_number(j) 
               + encode_number(k) + ")";
}

std::size_t ijk::shell_size() const {
    return ::shell_size(sum());
}

std::size_t ijk::index() const {
    return ijk_index(i, j, k);
}

std::size_t  ijk::operator[](std::size_t pos) const {
    assert(pos < 3);
    if (pos == 0) return i;
    if (pos == 1) return j;
    return k;
}

std::size_t& ijk::operator[](std::size_t pos) {
    assert(pos < 3);
    if (pos == 0) return i;
    if (pos == 1) return j;
    return k;
}


// relational operator
bool operator< (const ijk &a, const ijk &b) {
    const std::size_t aSum = a.sum();
    const std::size_t bSum = b.sum();

    if (aSum != bSum) return aSum < bSum;
    if (a.i != b.i) return a.i > b.i;
    if (a.j != b.j) return a.j > b.j;
    return a.k > b.k;
}

bool operator> (const ijk &a, const ijk &b) {
    return !(a < b) && !(a == b);
}

bool operator<=(const ijk &a, const ijk &b) {
    return (a < b) || (a == b);
}

bool operator>=(const ijk &a, const ijk &b) {
    return !(a < b);
}

bool operator==(const ijk &a, const ijk &b) {
    return (a.i == b.i) && (a.j == b.j) && (a.k == b.k);
}

bool operator!=(const ijk &a, const ijk &b) {
    return !(a == b);
}

std::vector<ijk> generate_shell(std::size_t n) {
    std::set<ijk> ijkSet;

    // TODO: new algorithm: listing ?
    for (std::size_t i = 0; i <= n; ++i) {
    for (std::size_t j = 0; j <= n; ++j) {
    for (std::size_t k = 0; k <= n; ++k) {
        if (i + j + k == n) {
            ijkSet.insert(ijk(i, j, k));
        }
    }}}

    return std::vector<ijk>(ijkSet.begin(), ijkSet.end());
}