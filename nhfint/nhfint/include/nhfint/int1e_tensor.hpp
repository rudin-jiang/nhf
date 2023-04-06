#pragma once

#include <vector>
#include <cstddef>

namespace nhfInt {

// plain one-electron integral of basis set
class Int1eTensorD0 {
public:
    std::size_t nBasis;
    std::size_t pageSize;

    std::vector<double> e1D0;   // nBasis * nBasis

    Int1eTensorD0();
    Int1eTensorD0(std::size_t nBasis);

    double  d0(std::size_t i, std::size_t j) const;
    double& d0(std::size_t i, std::size_t j);
};


class Int1eTensorD1 {
public:
    std::size_t nBasis;
    std::size_t pageSize;

    /*
     *  Data ordering: d0
     */
    std::vector<double> e1D1x;  // nBasis * nBasis
    std::vector<double> e1D1y;  // nBasis * nBasis
    std::vector<double> e1D1z;  // nBasis * nBasis

    Int1eTensorD1();
    Int1eTensorD1(std::size_t nBasis);

    // derivative of first basis function
    double  d1x0(std::size_t i, std::size_t j) const;
    double& d1x0(std::size_t i, std::size_t j);
    double  d1y0(std::size_t i, std::size_t j) const;
    double& d1y0(std::size_t i, std::size_t j);
    double  d1z0(std::size_t i, std::size_t j) const;
    double& d1z0(std::size_t i, std::size_t j);

    // specify the basis function to be derived
    double  d1x(std::size_t a, std::size_t i, std::size_t j) const;
    double& d1x(std::size_t a, std::size_t i, std::size_t j);
    double  d1y(std::size_t a, std::size_t i, std::size_t j) const;
    double& d1y(std::size_t a, std::size_t i, std::size_t j);
    double  d1z(std::size_t a, std::size_t i, std::size_t j) const;
    double& d1z(std::size_t a, std::size_t i, std::size_t j);
};


class Int1eTensorD2 {
public:
    std::size_t nBasis;
    std::size_t pageSize;

    /*
     *  Data ordering: d00 d01
     */
    std::vector<double> e1D2xx;     // 2 * nBasis * nBasis
    std::vector<double> e1D2yy;     // 2 * nBasis * nBasis
    std::vector<double> e1D2zz;     // 2 * nBasis * nBasis
    std::vector<double> e1D2xy;     // 2 * nBasis * nBasis
    std::vector<double> e1D2yz;     // 2 * nBasis * nBasis
    std::vector<double> e1D2zx;     // 2 * nBasis * nBasis

    Int1eTensorD2();
    Int1eTensorD2(std::size_t nBasis);

    double  d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j);
    double  d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j) const;
    double& d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j);

private:
    std::size_t _d2_index(
        std::size_t a, std::size_t b,
        std::size_t i, std::size_t j
    ) const;
};


} // namespace (nhfInt)



