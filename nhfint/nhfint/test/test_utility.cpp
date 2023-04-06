#include "nhfint/utility.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <vector>
#include <cstdio>

static const double absErr = 1e-3;

// Calculate the normalization constant using numerical integration
static double calc_gauss_norm(double a, std::size_t l) {
    double fwtm = 1.6651 / std::sqrt(a);

    // \int ( x^l * exp(-a * x^2) )^2  dx
    double xMax = 5.0 * fwtm;
    double dx   = 1e-4;
    double intVal = 0.0;
    for (double x = 0; x <= xMax; x += dx) {
        double gsVal = std::pow(x, l) * std::exp(-a*x*x);
        intVal += gsVal * gsVal;
    }

    return 1.0 / std::sqrt(2.0 * intVal * dx);
}

TEST(TestUtility, TestGaussNormCart1D) {
    for (double alpha = 0.1; alpha <= 5.0; alpha += 0.05) {
        for (std::size_t l = 0; l <= 8; ++l) {
            double norm0 = nhfInt::gauss_norm_cart(alpha, l);
            double norm1 = calc_gauss_norm(alpha, l);
            EXPECT_NEAR(norm0, norm1, absErr);
        }
    }
}

TEST(TestUtility, TestGaussNormCart3D) {
    for (double alpha = 0.1; alpha <= 5.0; alpha += 0.05) {
        for (std::size_t l = 0; l <= 8; ++l) {
        for (std::size_t m = 0; m <= 8; ++m) {
        for (std::size_t n = 0; n <= 8; ++n) {
            double norml = nhfInt::gauss_norm_cart(alpha, l);
            double normm = nhfInt::gauss_norm_cart(alpha, m);
            double normn = nhfInt::gauss_norm_cart(alpha, n);
            double norm0 = nhfInt::gauss_norm_cart(alpha, l, m , n);
            EXPECT_NEAR(norm0, norml * normm * normn, absErr);
        }}}
    }
}