#include "nhfint/calc_ericlass_dkr.hpp"
#include <gtest/gtest.h>
#include <unordered_map>
#include <unordered_set>
#include <cstddef>
#include <functional>
#include <random>
#include <vector>
#include <cmath>

// TODO:
// how to test float near?


static const double absErr = 1e-5;
static const double relErr = 1e-8;

static double relative_error(double val1, double val2) {
    return std::abs(val1 - val2) / std::abs(val2);
}

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


// random real generator
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(-10.0, 10.0);

static double random_real_gen() {
    return dis(gen);
}


// TODO:
// TEST(TestITensor, xxx) 


struct ITerm {
    std::size_t i, j, k, l;
};

bool operator==(const ITerm &a, const ITerm &b) {
    return (a.i == b.i) && (a.j == b.j)
            && (a.k == b.k) && (a.l == b.l);
}

struct iterm_hash {
    std::size_t operator()(const ITerm &a) const {
        return hash_val(a.i, a.j, a.k, a.l);
    }
};

using ItermValueMap = std::unordered_map<ITerm, double, iterm_hash>;
using ITermSet      = std::unordered_set<ITerm, iterm_hash>;


double dkr_hrr_bra_bfs(
    ItermValueMap &ivm, double x,
    std::size_t i, std::size_t j, std::size_t k, std::size_t l
) {
    auto it = ivm.find({i,j,k,l});

    if (it != ivm.end()) {
        return it -> second;
    }

    assert(j != 0);
    double ret = dkr_hrr_bra_bfs(ivm, x, i + 1, j - 1, k, l)
                    + x * dkr_hrr_bra_bfs(ivm, x, i, j - 1, k, l);
    
    // update ivm
    ivm[{i,j,k,l}] = ret;

    return ret;
}

double dkr_hrr_ket_bfs(
    ItermValueMap &ivm, double x,
    std::size_t i, std::size_t j, std::size_t k, std::size_t l
) {
    auto it = ivm.find({i,j,k,l});

    if (it != ivm.end()) {
        return it -> second;
    }

    assert(l != 0);
    double ret = dkr_hrr_ket_bfs(ivm, x, i, j, k + 1, l - 1)
                    + x * dkr_hrr_ket_bfs(ivm, x, i, j, k, l - 1);
    
    // update ivm
    ivm[{i,j,k,l}] = ret;

    return ret;
}

TEST(TestDkrHrr, TestDkrHrrBra) {
    std::vector<double> xList = {0.0, 0.5, 1.0, 5.0, 15.0};
    
    for (std::size_t na = 0; na <= 6; ++na) {
    for (std::size_t nb = 0; nb <= 6; ++nb) {
    for (std::size_t nc = 0; nc <= 6; ++nc) {
    for (std::size_t nd = 0; nd <= 6; ++nd) {
        std::size_t n = na + nb;
        nhfInt::dkr::ITensor n0xx(n, 0, nc, nd);

        for (std::size_t q = 0; q < n0xx.num0; ++q) {
            n0xx(q) = random_real_gen();
        }

        for (const double x : xList) {
            nhfInt::dkr::ITensor abxx = nhfInt::dkr::dkr_hrr_bra(n0xx, x, na, nb);

            ItermValueMap ivm;
            for (std::size_t i = 0; i < n0xx.cntA; ++i){
            for (std::size_t j = 0; j < n0xx.cntB; ++j){
            for (std::size_t k = 0; k < n0xx.cntC; ++k){
            for (std::size_t l = 0; l < n0xx.cntD; ++l){
                ivm[{i,j,k,l}] = n0xx(i,j,k,l);
            }}}}

            for (std::size_t i = 0; i < abxx.cntA; ++i){
            for (std::size_t j = 0; j < abxx.cntB; ++j){
            for (std::size_t k = 0; k < abxx.cntC; ++k){
            for (std::size_t l = 0; l < abxx.cntD; ++l){
                double bfsVal = dkr_hrr_bra_bfs(ivm, x, i, j, k, l);
                EXPECT_NEAR(abxx(i,j,k,l), bfsVal, absErr);
            }}}}
        }
    }}}}
}

TEST(TestDkrHrr, TestDkrHrrKet) {
    std::vector<double> xList = {0.0, 0.5, 1.0, 5.0, 15.0};
    
    for (std::size_t na = 0; na <= 6; ++na) {
    for (std::size_t nb = 0; nb <= 6; ++nb) {
    for (std::size_t nc = 0; nc <= 6; ++nc) {
    for (std::size_t nd = 0; nd <= 6; ++nd) {
        std::size_t m = nc + nd;
        nhfInt::dkr::ITensor xxm0(na, nb, m, 0);

        for (std::size_t q = 0; q < xxm0.num0; ++q) {
            xxm0(q) = random_real_gen();
        }

        for (const double x : xList) {
            nhfInt::dkr::ITensor xxcd = nhfInt::dkr::dkr_hrr_ket(xxm0, x, nc, nd);

            ItermValueMap ivm;
            for (std::size_t i = 0; i < xxm0.cntA; ++i){
            for (std::size_t j = 0; j < xxm0.cntB; ++j){
            for (std::size_t k = 0; k < xxm0.cntC; ++k){
            for (std::size_t l = 0; l < xxm0.cntD; ++l){
                ivm[{i,j,k,l}] = xxm0(i,j,k,l);
            }}}}

            for (std::size_t i = 0; i < xxcd.cntA; ++i){
            for (std::size_t j = 0; j < xxcd.cntB; ++j){
            for (std::size_t k = 0; k < xxcd.cntC; ++k){
            for (std::size_t l = 0; l < xxcd.cntD; ++l){
                double bfsVal = dkr_hrr_ket_bfs(ivm, x, i, j, k, l);
                EXPECT_NEAR(xxcd(i,j,k,l), bfsVal, absErr);
            }}}}
        }
    }}}}
}


double dkr_hrr_bfs(
    ItermValueMap &ivm, double xij, double xkl,
    std::size_t i, std::size_t j, std::size_t k, std::size_t l
) {
    auto it = ivm.find({i,j,k,l});

    if (it != ivm.end()) {
        return it -> second;
    }

    bool isFind = false;
    double ret = 0.0;

    if (j != 0 && isFind == false) {
        ret = dkr_hrr_bfs(ivm, xij, xkl, i + 1, j - 1, k, l)
                + xij * dkr_hrr_bfs(ivm, xij, xkl, i,     j - 1, k, l);
        isFind = true;
    }
    
    if (l != 0 && isFind == false){
        ret = dkr_hrr_bfs(ivm, xij, xkl, i, j, k + 1, l - 1)
                + xkl * dkr_hrr_bfs(ivm, xij, xkl, i, j, k,     l - 1);
        isFind = true;
    }

    assert(isFind);

    ivm[{i,j,k,l}] = ret;
    return ret;
}


TEST(TestDkrHrr, TestDkrHrr) {
    std::vector<double> xijList = {0.0, 0.5, 1.0, 5.0, 15.0};
    std::vector<double> xklList = {0.0, 0.5, 1.0, 5.0, 15.0};
    
    for (std::size_t na = 0; na <= 6; ++na) {
    for (std::size_t nb = 0; nb <= 6; ++nb) {
    for (std::size_t nc = 0; nc <= 6; ++nc) {
    for (std::size_t nd = 0; nd <= 6; ++nd) {
        std::size_t n = na + nb;
        std::size_t m = nc + nd;
        nhfInt::dkr::ITensor n0m0(n, 0, m, 0);

        for (std::size_t q = 0; q < n0m0.num0; ++q) {
            n0m0(q) = random_real_gen();
        }

        for (const double xij : xijList) {
        for (const double xkl : xklList) {

            nhfInt::dkr::ITensor abcd = nhfInt::dkr::dkr_hrr(n0m0, xij, xkl, na, nb, nc, nd);

            ItermValueMap ivm;
            for (std::size_t i = 0; i < n0m0.cntA; ++i){
            for (std::size_t j = 0; j < n0m0.cntB; ++j){
            for (std::size_t k = 0; k < n0m0.cntC; ++k){
            for (std::size_t l = 0; l < n0m0.cntD; ++l){
                ivm[{i,j,k,l}] = n0m0(i,j,k,l);
            }}}}

            for (std::size_t i = 0; i < abcd.cntA; ++i){
            for (std::size_t j = 0; j < abcd.cntB; ++j){
            for (std::size_t k = 0; k < abcd.cntC; ++k){
            for (std::size_t l = 0; l < abcd.cntD; ++l){
                double bfsVal = dkr_hrr_bfs(ivm, xij, xkl, i, j, k, l);
                double funVal = abcd(i,j,k,l);
                if (std::abs(bfsVal) < absErr || std::abs(funVal) < absErr) {
                    EXPECT_NEAR(funVal, bfsVal, absErr);
                }
                else {
                    EXPECT_TRUE(relative_error(bfsVal, funVal) < relErr);
                }
            }}}}
        }}
    }}}}
}



