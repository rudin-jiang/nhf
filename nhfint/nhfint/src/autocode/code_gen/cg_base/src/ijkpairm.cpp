#include "cg_base/ijkpairm.hpp"
#include "cg_base/ijk.hpp"
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <cassert>

ijkPairM::ijkPairM(): m(0) {}

ijkPairM::ijkPairM(ijk a, ijk b, std::size_t m)
: a(a), b(b), m(m) {}


std::string ijkPairM::to_string() const {
    char buff[128];
    std::sprintf(buff, "%2zu  %2zu  %2zu  %2zu  %2zu  %2zu  %2zu",
                                 a.i, a.j, a.k, b.i, b.j, b.k, m);
    return std::string(buff);
}

std::size_t ijkPairM::single_pair_size() const {
    return a.shell_size() * b.shell_size();
}

std::size_t ijkPairM::single_index() const {
    return a.index() * b.shell_size() + b.index();
}

std::size_t ijkPairM::cnt_zero() const {
    return a.cnt_zero() + b.cnt_zero();
}

std::size_t ijkPairM::cnt_nonz() const {
    return a.cnt_nonz() + b.cnt_nonz();
}

std::size_t  ijkPairM::operator[](std::size_t pos) const {
    assert(pos < 6);
    if (pos < 3) return a[pos];
    return b[pos - 3];
}

std::size_t& ijkPairM::operator[](std::size_t pos) {
    assert(pos < 6);
    if (pos < 3) return a[pos];
    return b[pos - 3];
}

std::vector<std::size_t> ijkPairM::nonzero_index() const {
    std::vector<std::size_t> ret;
    if (a.i != 0) ret.push_back(0);
    if (a.j != 0) ret.push_back(1);
    if (a.k != 0) ret.push_back(2);
    if (b.i != 0) ret.push_back(3);
    if (b.j != 0) ret.push_back(4);
    if (b.k != 0) ret.push_back(5);
    return ret;
}

std::vector<ijkPairM> ijkPairM::hgp_vrr_need(std::size_t pos) const {
    assert((*this)[pos] != 0);
    
}


bool operator< (const ijkPairM &ab, const ijkPairM &cd) {
    std::size_t abSumA = ab.a.sum();
    std::size_t abSumB = ab.b.sum();
    std::size_t cdSumA = cd.a.sum();
    std::size_t cdSumB = cd.b.sum();

    if (abSumA != cdSumA) return abSumA < cdSumA;
    if (abSumB != cdSumB) return abSumB < cdSumB;

    if (ab.m != cd.m) return ab.m < cd.m;

    if (ab.a != cd.a) return ab.a < cd.a;
    return ab.b < cd.b;
}

bool operator> (const ijkPairM &ab, const ijkPairM &cd) {
    return !(ab == cd) && !(ab < cd);
}

bool operator<=(const ijkPairM &ab, const ijkPairM &cd) {
    return (ab == cd) || (ab < cd);
}

bool operator>=(const ijkPairM &ab, const ijkPairM &cd) {
    return !(ab < cd);
}

bool operator==(const ijkPairM &ab, const ijkPairM &cd) {
    return (ab.a == cd.a) && (ab.b == cd .b) && (ab.m == cd.m);
}

bool operator!=(const ijkPairM &ab, const ijkPairM &cd) {
    return !(ab == cd);
}

std::vector<ijkPairM> generate_ijkpair0(std::size_t na, std::size_t nb) {
    std::vector<ijk> ijkVecA = generate_shell(na);
    std::vector<ijk> ijkVecB = generate_shell(nb);

    std::set<ijkPairM> pairmSet;
    for (const ijk &a : ijkVecA) {
    for (const ijk &b : ijkVecB) {
        pairmSet.insert(ijkPairM(a, b, 0));
    }}

    return std::vector<ijkPairM> (
        pairmSet.begin(), pairmSet.end()
    );
}
