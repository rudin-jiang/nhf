#include "cg_base/iterm.hpp"
#include <cstddef>
#include <vector>
#include <set>
#include <string>
#include <cstdio>
#include <cassert>


ITerm::ITerm(): i(0), j(0), k(0), l(0) {}

ITerm::ITerm(std::size_t i, std::size_t j,
             std::size_t k, std::size_t l)
: i(i), j(j), k(k), l(l) {
    assert(i < 100);
    assert(j < 100);
    assert(k < 100);
    assert(l < 100);
}

std::string ITerm::to_string() const {
    char buff[64];
    std::sprintf(buff, "%2zu  %2zu  %2zu  %2zu",
                                        i, j, k, l);
    return std::string(buff);
}


bool operator<(const ITerm &a, const ITerm &b) {
    if (a.i != b.i) return a.i < b.i;
    if (a.j != b.j) return a.j < b.j;
    if (a.k != b.k) return a.k < b.k;
    return a.l < b.l;
}

bool operator>(const ITerm &a, const ITerm &b) {
    return !(a < b) && !(a == b);
}

bool operator<=(const ITerm &a, const ITerm &b) {
    return (a < b) || (a == b);
}

bool operator>=(const ITerm &a, const ITerm &b) {
    return !(a < b);
}

bool operator==(const ITerm &a, const ITerm &b) {
    return (a.i == b.i) && (a.j == b.j)
           && (a.k == b.k) && (a.l == b.l);
}

bool operator!=(const ITerm &a, const ITerm &b) {
    return !(a == b);
}


std::vector<ITerm> generate_iterm(std::size_t na, std::size_t nb, 
                                  std::size_t nc, std::size_t nd) {
    assert(na < 100);
    assert(nb < 100);
    assert(nc < 100);
    assert(nd < 100);

    const std::size_t nITerm = (na+1)*(nb+1)*(nc+1)*(nd+1);

    std::vector<ITerm> iTermVec;
    iTermVec.reserve(nITerm);

    std::size_t n = 0;
    for (std::size_t i = 0; i <= na; ++i) {
    for (std::size_t j = 0; j <= nb; ++j) {
    for (std::size_t k = 0; k <= nc; ++k) {
    for (std::size_t l = 0; l <= nd; ++l) {
        assert(n < nITerm);
        iTermVec.push_back(ITerm(i, j, k, l));
    }}}}

    return iTermVec;
}