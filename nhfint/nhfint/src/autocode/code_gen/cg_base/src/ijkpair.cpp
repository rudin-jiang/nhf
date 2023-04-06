#include "cg_base/ijkpair.hpp"
#include "cg_base/ijk.hpp"
#include "cg_base/utility.hpp"
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

std::string ijkPair::to_string() const {
    char buff[64];
    std::sprintf(buff, "%2zu  %2zu  %2zu  %2zu  %2zu %2zu",
                        a.i, a.j, a.k, b.i, b.j, b.k);
    return std::string(buff);
}

std::string ijkPair::to_comment() const {
    return "[" + a.to_comment() + b.to_comment() + "]";
}

std::size_t ijkPair::pair_size() const {
    return a.shell_size() * b.shell_size();
}

std::size_t ijkPair::index() const {
    return a.index() * b.shell_size() + b.index();
}

// form in axbx[idx]
std::string ijkPair::array_element(std::size_t indexLen) const {
    return array_name(a.sum(), b.sum()) + "[" + format_number(index(), indexLen) + "]";
}

std::size_t  ijkPair::operator[](std::size_t pos) const {
    assert(pos < 6);
    if (pos < 3) return a[pos];
    return b[pos - 3];
}

std::vector<ijkPair> ijkPair::hgp_hrr_need(std::size_t pos) const {
    assert(b[pos] != 0);

    ijkPair needA = (*this);
    ijkPair needB = (*this);

    needA.a[pos] += 1;
    needA.b[pos] -= 1;
    needB.b[pos] -= 1;

    return {needA, needB};
}

std::size_t& ijkPair::operator[](std::size_t pos) {
    assert(pos < 6);
    if (pos < 3) return a[pos];
    return b[pos - 3];
}

bool operator< (const ijkPair &ab, const ijkPair &cd) {
    std::size_t abASum = ab.a.sum();
    std::size_t abBSum = ab.b.sum();
    std::size_t cdASum = cd.a.sum();
    std::size_t cdBSum = cd.b.sum();

    if (abASum != cdASum) return abASum < cdASum;
    if (abBSum != cdBSum) return abBSum < cdBSum;

    if (ab.a != cd.a) return ab.a < cd.a;
    return ab.b < cd.b;
}

bool operator> (const ijkPair &ab, const ijkPair &cd) {
    return !(ab < cd) && !(ab == cd);
}

bool operator<=(const ijkPair &ab, const ijkPair &cd) {
    return (ab < cd) || (ab == cd);
}

bool operator>=(const ijkPair &ab, const ijkPair &cd) {
    return !(ab < cd);
}

bool operator==(const ijkPair &ab, const ijkPair &cd) {
    return (ab.a == cd.a) && (ab.b == cd.b);
}

bool operator!=(const ijkPair &ab, const ijkPair &cd) {
    return !(ab == cd);
}

std::vector<ijkPair> generate_pair(std::size_t na, std::size_t nb) {
    std::set<ijkPair> pairSet;

    std::vector<ijk> ijkVecA = generate_shell(na);
    std::vector<ijk> ijkVecB = generate_shell(nb);

    // for (const ijk &a : ijkVecA) {
    // for (const ijk &b : ijkVecB) {
    //     pairSet.insert(ijkPair(a,b));
    // }}

    for (std::size_t i = 0; i < ijkVecA.size(); ++i) {
    for (std::size_t j = 0; j < ijkVecB.size(); ++j) {
        pairSet.insert(ijkPair(ijkVecA[i], ijkVecB[j]));
    }}


    return std::vector<ijkPair>(
        pairSet.begin(), pairSet.end()
    );
}
