#include "nhfint/calc_ericlass_tho.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <vector>
#include <cstdio>

const double eps = 1e-5;
const double relErr = 1e-5;

static double relative_error(double val0, double val1) {
    assert(std::abs(val0) > eps);
    assert(std::abs(val1) > eps);
    return std::abs(val0 - val1) / std::abs(val0);
}


TEST(TestUtility, TestBilomialPrefactor) {

    const double absErrTest = 1e-5;

    std:std::vector<double> xList = {-1.0, -0.5, 0.0, 0.5, 1.0};

    for (double a = -5.0; a < 5.0; a += 0.1) {
    for (double b = -5.0; b < 5.0; b += 0.1) {
        for (std::size_t l = 0; l <= 8; ++l){
        for (std::size_t m = 0; m <= 8; ++m){
            for (const double x : xList){
                double val0 = std::pow(x+a, l) * std::pow(x+b, m);
                double val1 = 0.0;
                for (std::size_t j = 0; j <= l+m; ++j) {
                    val1 += nhfInt::tho::binomial_prefactor(j, l, m, a,b)
                            * std::pow(x, j);
                }
                
                if (std::abs(val0) <= eps || std::abs(val1) <= eps) {
                    EXPECT_NEAR(val0, val1, 0.1*eps);
                }
                else {
                    EXPECT_TRUE(relative_error(val0, val1) < relErr);
                }

                if (relative_error(val0, val1) >= relErr)
                std::printf("a = %.2lf    b = %.2lf    l = %zu    m = %zu   x = %.2lf \n", a, b, l, m, x);
            }
        }}
    }}
}
