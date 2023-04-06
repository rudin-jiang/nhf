#include "cg_base/ijk.hpp"
#include <gtest/gtest.h>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>
#include <random>

TEST(TestIjk, TestConstructors) {
    // default
    ijk a0;
    EXPECT_EQ(a0.i, 0);
    EXPECT_EQ(a0.j, 0);
    EXPECT_EQ(a0.k, 0);

    // value
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a1(i, j, k);
        EXPECT_EQ(a1.i, i);
        EXPECT_EQ(a1.j, j);
        EXPECT_EQ(a1.k, k);
    }}}

    // string
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        char buff[64];
        std::sprintf(buff, "%zu %zu %zu", i, j, k);
        
        std::string line(buff);
        ijk a2(line);
        EXPECT_EQ(a2.i, i);
        EXPECT_EQ(a2.j, j);
        EXPECT_EQ(a2.k, k);
    }}}
}

TEST(TestIjk, TestSum) {
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a(i, j, k);
        EXPECT_EQ(a.sum(), i + j + k);
    }}}
}

TEST(TestIjk, TestCntZero) {
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a(i, j, k);

        std::size_t nZero = 0;
        if (i == 0) ++nZero;
        if (j == 0) ++nZero;
        if (k == 0) ++nZero;

        EXPECT_EQ(a.cnt_zero(), nZero);
    }}}
}

TEST(TestIjk, TestCntNonz) {
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a(i, j, k);

        std::size_t nNonZero = 0;
        if (i != 0) ++nNonZero;
        if (j != 0) ++nNonZero;
        if (k != 0) ++nNonZero;

        EXPECT_EQ(a.cnt_nonz(), nNonZero);
    }}}
}

TEST(TestIjk, TestNonzeroIndex) {
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        std::vector<std::size_t> nonzIdx;
        if (i != 0) nonzIdx.push_back(0);
        if (j != 0) nonzIdx.push_back(1);
        if (k != 0) nonzIdx.push_back(2);

        ijk a(i, j, k);
        EXPECT_EQ(a.nonzero_index(), nonzIdx);
    }}}
}

TEST(TestIjk, TestToString) {
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a(i, j, k);
        std::string iStr = (i < 10 ? " " : "") + std::to_string(i);
        std::string jStr = (j < 10 ? " " : "") + std::to_string(j);
        std::string kStr = (k < 10 ? " " : "") + std::to_string(k);
        std::string line = iStr + "  " + jStr + "  " + kStr;
        EXPECT_EQ(a.to_string(), line);
    }}}
}

TEST(TestIjk, TestIndex) {
    for (std::size_t n = 0; n <= 3; ++n) {
        std::vector<ijk> ijkVec = generate_shell(n);
        std::size_t ijkNumCalc = (n+1)*(n+2)/2;
        EXPECT_EQ(ijkVec.size(), ijkNumCalc);

        for (std::size_t i = 0; i < ijkVec.size(); ++i) {
            EXPECT_EQ(ijkVec[i].index(), i);
        }
    }
}

TEST(TestIjk, TestOperatorSquareBracket) {
    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a(i, j, k);
        EXPECT_EQ(a[0], i);
        EXPECT_EQ(a[1], j);
        EXPECT_EQ(a[2], k);
    }}}

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> distrib(0, 10);

    std::size_t newI = distrib(gen);
    std::size_t newJ = distrib(gen);
    std::size_t newK = distrib(gen);

    for (std::size_t i = 0; i <= 10; ++i) {
    for (std::size_t j = 0; j <= 10; ++j) {
    for (std::size_t k = 0; k <= 10; ++k) {
        ijk a(i, j, k);

        a[0] = newI;
        EXPECT_EQ(a[0], newI);
        EXPECT_EQ(a[1], j);
        EXPECT_EQ(a[2], k);

        a[1] = newJ;
        EXPECT_EQ(a[0], newI);
        EXPECT_EQ(a[1], newJ);
        EXPECT_EQ(a[2], k);

        a[2] = newK;
        EXPECT_EQ(a[0], newI);
        EXPECT_EQ(a[1], newJ);
        EXPECT_EQ(a[2], newK);
    }}}
}
