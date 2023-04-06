#include "boysfun/boysfun.hpp"
#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>

// When the real value is less than eps,
// the value can be considered as 0.
static const double eps = 1e-99;

// When the actual value is greater than eps,
// the relative error between the calculated 
// value and the actual value is required to
// be less than this value (relErr).
static const double relErr = 1e-13;
static const double absErr = 1e-15;


static double relative_error(double calcVal, double realVal) {
    assert(realVal >= eps);
    return std::abs(realVal - calcVal) / std::abs(realVal);
}


static void test_boysfun(const std::string &filename) {
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;

    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
    }

    for (const double x : xValueList) {
        for (std::size_t n = 0; n <= nhfInt::boysFun::BoysFunMaxN; ++n) {
            // read in real value
            double realVal = 0.0;
            ifs >> realVal;
            assert(realVal >= 0.0);

            // calculated value
            double calcVal = nhfInt::boysFun::boysfun(n, x);

            // boysfun value is always positive
            EXPECT_TRUE(calcVal >= 0.0);

            // relative error test
            if (realVal < eps) {
                EXPECT_TRUE(calcVal < realVal * 1.01);
            }
            else {
                double relErrNow = relative_error(calcVal, realVal);
                EXPECT_TRUE(relErrNow < relErr);
            }

            // absolute error test
            EXPECT_NEAR(calcVal, realVal, absErr);
        }
    }
}

TEST(TestBoysfun, TestBoysfunSmallX) {
    std::string filename = "test_data/boysfun_test_data_small_x.dat";
    test_boysfun(filename);
}

TEST(TestBoysfun, TestBoysfunLargeX) {
    std::string filename = "test_data/boysfun_test_data_large_x.dat";
    test_boysfun(filename);
}