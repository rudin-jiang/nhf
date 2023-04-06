#include <cstddef>
#define _USE_MATH_DEFINES

#include "nhfint/utility.hpp"
#include "nhfmath/mathfun.hpp"
#include <cmath>
#include <cstdint>



namespace nhfInt {

double gauss_norm_cart(double alpha, std::size_t l) {
    double lSF = nhfMath::semifactorial(2 * std::int64_t(l) - 1);
    return std::pow(M_2_PI * alpha, 0.25)
            * std::sqrt(std::pow(4.0 * alpha, l) / lSF);
}

double gauss_norm_cart(double alpha, std::size_t l, std::size_t m, std::size_t n) {
    double lSF = nhfMath::semifactorial(2 * std::int64_t(l) - 1);
    double mSF = nhfMath::semifactorial(2 * std::int64_t(m) - 1);
    double nSF = nhfMath::semifactorial(2 * std::int64_t(n) - 1);
    return std::pow(M_2_PI * alpha, 0.75)
            * std::sqrt(std::pow(4.0 * alpha, l+m+n) / (lSF*mSF*nSF));
}




} // namespace (nhfInt)

