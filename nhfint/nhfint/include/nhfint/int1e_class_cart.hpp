#pragma once

#include <cstddef>
#include <vector>


namespace nhfInt {

class Int1eClassCartD0 {
public:
    std::size_t ang0, ang1;     // angular moment of basis
    std::size_t nbs0, nbs1;     // (l+1) * (l+2) / 2
    std::size_t num0, num1;     // for fast indexing
    std::vector<double> e1D0;   // row-major matrix

    Int1eClassCartD0();
    Int1eClassCartD0(std::size_t ang0, std::size_t ang1);

    double  d0(std::size_t i, std::size_t j) const;
    double& d0(std::size_t i, std::size_t j);
};


class Int1eClassCartD1 {
public:
    std::size_t ang0, ang1;     // angular moment of basis
    std::size_t nbs0, nbs1;     // (l+1) * (l+2) / 2
    std::size_t num0, num1;     // for fast indexing
    
    /*
     * Derivatives of first basis function are
     * placed in the first half of the vectors.
     *  - d0 d1
     */
    std::vector<double> e1D1x;    // size = 2 * D0
    std::vector<double> e1D1y;    // size = 2 * D0
    std::vector<double> e1D1z;    // size = 2 * D0

    Int1eClassCartD1();
    Int1eClassCartD1(std::size_t ang0, std::size_t ang1);

    // Derivatives of first basis function
    double  d1x0(std::size_t i, std::size_t j) const;
    double& d1x0(std::size_t i, std::size_t j);
    double  d1y0(std::size_t i, std::size_t j) const;
    double& d1y0(std::size_t i, std::size_t j);
    double  d1z0(std::size_t i, std::size_t j) const;
    double& d1z0(std::size_t i, std::size_t j);

    // Derivatives of second basis function
    double  d1x1(std::size_t i, std::size_t j) const;
    double& d1x1(std::size_t i, std::size_t j);
    double  d1y1(std::size_t i, std::size_t j) const;
    double& d1y1(std::size_t i, std::size_t j);
    double  d1z1(std::size_t i, std::size_t j) const;
    double& d1z1(std::size_t i, std::size_t j);

    // Specify the basis function to be derived
    double  d1x(std::size_t a, std::size_t i, std::size_t j) const;
    double& d1x(std::size_t a, std::size_t i, std::size_t j);
    double  d1y(std::size_t a, std::size_t i, std::size_t j) const;
    double& d1y(std::size_t a, std::size_t i, std::size_t j);
    double  d1z(std::size_t a, std::size_t i, std::size_t j) const;
    double& d1z(std::size_t a, std::size_t i, std::size_t j);
};


class Int1eClassCartD2 {
public:
    std::size_t ang0, ang1;     // angular moment of basis
    std::size_t nbs0, nbs1;     // (l+1) * (l+2) / 2
    std::size_t num0, num1;     // for fast indexing
    
    /*
     * Data ordering in vectors: 
     *  - d00 d01 d10 d11
     */
    std::vector<double> e1D2xx;   // size = 2 * 2 * D0
    std::vector<double> e1D2yy;   // size = 2 * 2 * D0
    std::vector<double> e1D2zz;   // size = 2 * 2 * D0
    std::vector<double> e1D2xy;   // size = 2 * 2 * D0
    std::vector<double> e1D2yz;   // size = 2 * 2 * D0
    std::vector<double> e1D2zx;   // size = 2 * 2 * D0

    Int1eClassCartD2();
    Int1eClassCartD2(std::size_t ang0, std::size_t ang1);

    // Specify the basis function to be derived
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
};


} // namespace (nhfInt)

