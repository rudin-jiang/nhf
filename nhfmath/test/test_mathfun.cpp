#include "nhfmath/mathfun.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <cstddef>
#include <cstdint>

static const double absErr = 1.0e-14;

TEST(TestMathFun, TestFactorial) {
    for (std::size_t n = 0; n <= nhfMath::FactorialMaxN; ++n) {
        double factorialCalc = 1.0;
        for (std::size_t i = 1; i <= n; ++i) {
            factorialCalc *= i;
        }

        EXPECT_NEAR(nhfMath::factorial(n), factorialCalc, absErr);
    }
}


TEST(TestMathFun, TestCombination) {
    for (std::size_t n = 0; n <= nhfMath::CombinationMaxN; ++n) {
        for(std::size_t k = 0; k <= n; ++k) {

            double combinationCalc = 1.0;
            for (std::size_t t = n, i = 1; i <= k; ++i, --t) {
                combinationCalc *= t;
                combinationCalc /= i;
            }

            EXPECT_NEAR(nhfMath::combination(n, k), combinationCalc, absErr);
        }
    }
}


TEST(TestMathFun, TestSemifactorial) {
    for (std::int64_t n = nhfMath::SemifactorialMinN;
                        n <= nhfMath::SemifactorialMaxN; ++n) {
        double semifactorialCalc = 1.0;
        if (n >= 0) {
            for (std::int64_t i = ((n&1) ? 1 : 2); i <= n; i += 2){
                semifactorialCalc *= i;
            }
        }

        EXPECT_NEAR(nhfMath::semifactorial(n), semifactorialCalc, absErr);
    }
}