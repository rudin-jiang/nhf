// 2023-01-15
// TODO: need to be reviewed
//       translate comments

/*  
 *  1、为了复用中间变量，总是同时计算一个class
 *  2、使用MD方法计算双电子积分时，需要知道组合系数Eijt和辅助函数Rijk
 *  3、这个小程序验证了在计算同一个双电子积分class的时候，需要哪些Eijt和Rijk
 *  4、结论：给定一个ericlass，四个角动量分别为 na，nb，nc，nd，
 *          Eijt: i = 0 ~ na+nb, j = 0 ~ nc+nd, t = 0 ~ i+j
 *          Rijk: 0 <= i+j+k <= na+nb+nc+nd
 */

// #define NDEBUG
#include <iostream>
#include <cassert>
#include <vector>
#include <cstddef>
#include <functional>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <cstdio>


// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html
template <typename T>
inline void hash_combine(std::size_t &seed, const T &val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// auxiliary generic functions to create a hash value using a seed
template <typename T>
inline void hash_val(std::size_t &seed, const T &val) {
    hash_combine(seed, val);
}

template <typename T, typename... Types>
inline void hash_val(std::size_t &seed, const T &val, const Types &... args) {
    hash_combine(seed, val);
    hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}


// class ijk:
// 1. record angular momentum of a basis function
// 2. record auxiliary generic functions R(N,M,L) needed in MD scheme
// see ref: doi: 10.1016/0021-9991(78)90092-X
struct ijk {
    std::size_t i, j, k;

    ijk() : i(0), j(0), k(0) {}

    ijk(std::size_t i, std::size_t j, std::size_t k)
    : i(i), j(j), k(k) {}

    std::size_t sum() const { return i + j + k; }
};

bool operator< (const ijk &a, const ijk &b);
bool operator> (const ijk &a, const ijk &b);
bool operator<=(const ijk &a, const ijk &b);
bool operator>=(const ijk &a, const ijk &b);
bool operator==(const ijk &a, const ijk &b);
bool operator!=(const ijk &a, const ijk &b);

// Logic of relational operators
bool operator<(const ijk &a, const ijk &b) {
    const std::size_t aSum = a.sum();
    const std::size_t bSum = b.sum();
    if (aSum != bSum) return aSum < bSum;

    if (a.i != b.i) return a.i < b.i;
    if (a.j != b.j) return a.j < b.j;
    return a.k < b.k;
}

bool operator>(const ijk &a, const ijk &b) {
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

// hash function of struct ijk
struct ijk_hash {
    std::size_t operator()(const ijk &a) const {
        return hash_val(a.i, a.j, a.k);
    }
};


// generate all ijk that i+j+k == n
// and i >= 0 && j >= 0 && k >= 0
// TODO: more efficient algorithm
std::vector<ijk> generate_ijk(std::size_t n) {
    std::set<ijk> ijkSet;
    for (std::size_t i = 0; i <= n; ++i) {
    for (std::size_t j = 0; j <= n; ++j) {
    for (std::size_t k = 0; k <= n; ++k) {
        if (i + j + k == n) {
            ijkSet.insert(ijk(i, j, k));
        }
    }}}

    // check the number of all Rijk
    const std::size_t numCalc = (n+1) * (n+2) / 2;
    assert(ijkSet.size() == numCalc);

    return std::vector<ijk>(ijkSet.begin(), ijkSet.end());
}


// class eii:
// 
struct eii {
    std::size_t i, j, t;

    eii() : i(0), j(0), t(0) {}

    eii(std::size_t i, std::size_t j, std::size_t t)
    : i(i), j(j), t(t) {}
};

bool operator< (const eii &a, const eii &b);
bool operator> (const eii &a, const eii &b);
bool operator<=(const eii &a, const eii &b);
bool operator>=(const eii &a, const eii &b);
bool operator==(const eii &a, const eii &b);
bool operator!=(const eii &a, const eii &b);

// Logic of relational operators
bool operator<(const eii &a, const eii &b) {
    const std::size_t aSum = a.i + a.j;
    const std::size_t bSum = b.i + b.j;
    if (aSum != bSum) return aSum < bSum;
    if (a.i != b.i) return a.i < b.i;
    return a.t < b.t;
}

bool operator>(const eii &a, const eii &b) {
    return !(a < b) && !(a == b);
}

bool operator<=(const eii &a, const eii &b) {
    return (a < b) || (a == b);
}

bool operator>=(const eii &a, const eii &b) {
    return ! (a < b);
}

bool operator==(const eii &a, const eii &b) {
    return (a.i == b.i) && (a.j == b.j) && (a.t == b.t);
}

bool operator!=(const eii &a, const eii &b) {
    return ! (a == b);
}

// hash function of struct eii
struct eii_hash {
    std::size_t operator()(const eii &a) const {
        return hash_val(a.i, a.j, a.t);
    }
};

// generate all eii that 0 <= t <= i+j
std::vector<eii> generate_eii(std::size_t i, std::size_t j) {
    std::vector<eii> eiiVec(i+j+1);
    for (std::size_t t = 0; t <= i+j; ++t) {
        eiiVec[t] = eii(i, j, t);
    }

    return eiiVec;
}

// given na, nb, nc, nd
// return all Eijt and Rijk needed
// using eii record Eijt
// using ijk record Rijk
void eijt_rijk_need(
    std::vector<eii> &eijtNeed,
    std::vector<ijk> &rijkNeed,
    std::size_t na, std::size_t nb,
    std::size_t nc, std::size_t nd
) {
    std::unordered_set<eii, eii_hash>   eijtSet;
    std::unordered_set<ijk, ijk_hash>   rijkSet;

    std::vector<ijk> bsSetA = generate_ijk(na);
    std::vector<ijk> bsSetB = generate_ijk(nb);
    std::vector<ijk> bsSetC = generate_ijk(nc);
    std::vector<ijk> bsSetD = generate_ijk(nd);

    for (const ijk &a : bsSetA) {
    for (const ijk &b : bsSetB) {
    for (const ijk &c : bsSetC) {
    for (const ijk &d : bsSetD) {
        for (std::size_t ta = 0; ta <= a.i + b.i; ++ta) {
        for (std::size_t ua = 0; ua <= a.j + b.j; ++ua) {
        for (std::size_t va = 0; va <= a.k + b.k; ++va) {
        for (std::size_t tb = 0; tb <= c.i + d.i; ++tb) {
        for (std::size_t ub = 0; ub <= c.j + d.j; ++ub) {
        for (std::size_t vb = 0; vb <= c.k + d.k; ++vb) {
            eijtSet.insert(eii(a.i, b.i, ta));
            rijkSet.insert(ijk(ta+tb, ua+ub, va+vb));
        }}}}}}
    }}}}

    eijtNeed = std::vector<eii>(eijtSet.begin(), eijtSet.end());
    rijkNeed = std::vector<ijk>(rijkSet.begin(), rijkSet.end());

    std::sort(eijtNeed.begin(), eijtNeed.end());
    std::sort(rijkNeed.begin(), rijkNeed.end());
}


// guess eijt need
// i range from 0 to na
// j range from 0 to nb
// t range from 0 to i + j
std::vector<eii> guess_eijt_need(std::size_t na, std::size_t nb) {
    std::vector<eii> ret;
    for (std::size_t a = 0; a <= na; ++a) {
    for (std::size_t b = 0; b <= nb; ++b) {
        std::vector<eii> eijtVec = generate_eii(a, b);
        for (const eii &e : eijtVec) {
            ret.push_back(e);
        }
    }}
    std::sort(ret.begin(), ret.end());
    return ret;
}


// guess rijk need
// all rijk that 0 <= i+j+k <= na+nb+nc+nd
std::vector<ijk> guess_rijk_need(std::size_t nSum) {
    std::vector<ijk> ret;
    for (std::size_t n = 0; n <= nSum; ++n) {
        std::vector<ijk> rijkVec = generate_ijk(n);
        for (const ijk &i : rijkVec) {
            ret.push_back(i);
        }
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}


int main() {
    const std::size_t nMax = 6;

    std::vector<eii> eijtNeed;
    std::vector<ijk> rijkNeed;

    for (std::size_t na = 0; na <= nMax; ++na) {
    for (std::size_t nb = 0; nb <= nMax; ++nb) {
    for (std::size_t nc = 0; nc <= nMax; ++nc) {
    for (std::size_t nd = 0; nd <= nMax; ++nd) {

        std::printf("na = %zu    nb = %zu    "
                    "nc = %zu    nd = %zu    ",
                    na, nb, nc, nd);

        eijt_rijk_need(eijtNeed, rijkNeed, na, nb, nc, nd);

        std::vector<eii> eijtGuess = guess_eijt_need(na, nb);
        std::vector<ijk> rijkGuess = guess_rijk_need(na + nb + nc + nd);

        if (eijtGuess == eijtNeed) {
            std::printf("eijt  ok       ");
        }
        else {
            std::printf("eijt  error    ");
        }

        if (rijkGuess == rijkNeed) {
            std::printf("rijk  ok       ");
        }
        else {
            std::printf("rijk  error    ");
        }

        std::cout << std::endl;
    }}}}

    return 0;
}