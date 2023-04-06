#pragma once

#include <cstddef>
#include <vector>

namespace nhfInt {

class Int2eClassCartD0 {
public:
    std::size_t ang0, ang1, ang2, ang3;     // angular moment of basis
    std::size_t nbs0, nbs1, nbs2, nbs3;     // (l+1) * (l+2) / 2
    std::size_t num0, num1, num2, num3;     // for fast indexing
    std::vector<double> e2D0;

    Int2eClassCartD0();
    Int2eClassCartD0(std::size_t ang0, std::size_t ang1,
                     std::size_t ang2, std::size_t ang3);

    double  d0(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d0(std::size_t i, std::size_t j, std::size_t k, std::size_t l);

    double  d0(std::size_t ij, std::size_t kl) const;
    double& d0(std::size_t ij, std::size_t kl);

    double  d0(std::size_t ijkl) const;
    double& d0(std::size_t ijkl);
};


class Int2eClassCartD1 {
public:
    std::size_t ang0, ang1, ang2, ang3;     // angular moment of basis
    std::size_t nbs0, nbs1, nbs2, nbs3;     // (l+1) * (l+2) / 2
    std::size_t num0, num1, num2, num3;     // for fast indexing

    /* TODO: better comments
     * Data ordering:
     *  - d0 d1 d2 d3
     */
    std::vector<double> e2D1x;    // size = 4 * D0
    std::vector<double> e2D1y;    // size = 4 * D0
    std::vector<double> e2D1z;    // size = 4 * D0

    Int2eClassCartD1();
    Int2eClassCartD1(std::size_t ang0, std::size_t ang1,
                     std::size_t ang2, std::size_t ang3);

    double  d1x0(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1x0(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1y0(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1y0(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1z0(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1z0(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);

    double  d1x1(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1x1(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1y1(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1y1(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1z1(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1z1(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);

    double  d1x2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1x2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1y2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1y2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1z2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1z2(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);

    double  d1x3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1x3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1y3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1y3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1z3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1z3(std::size_t i, std::size_t j, std::size_t k ,std::size_t l);

    // Specify the basis function to be derived
    double  d1x(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1x(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1y(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1y(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1z(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1z(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
};


class Int2eClassCartD2 {
public:
    std::size_t ang0, ang1, ang2, ang3;     // angular moment of basis
    std::size_t nbs0, nbs1, nbs2, nbs3;     // (l+1) * (l+2) / 2
    std::size_t num0, num1, num2, num3;     // for fast indexing

    /* TODO: better comments
     * Data ordering:
     *  - d00 d01 d02 d03
     *    d10 d11 d12 d13
     *    d20 d21 d22 d23
     *    d30 d31 d32 d33
     */
    std::vector<double> e2D2xx;   // size = 4 * 4 * D0
    std::vector<double> e2D2yy;   // size = 4 * 4 * D0
    std::vector<double> e2D2zz;   // size = 4 * 4 * D0
    std::vector<double> e2D2xy;   // size = 4 * 4 * D0
    std::vector<double> e2D2yz;   // size = 4 * 4 * D0
    std::vector<double> e2D2zx;   // size = 4 * 4 * D0

    Int2eClassCartD2();
    Int2eClassCartD2(std::size_t ang0, std::size_t ang1,
                     std::size_t ang2, std::size_t ang3);

    // Specify the basis function to be derived
    double  d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2xx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2yy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2zz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2xy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2yz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2zx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2yx(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2zy(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d2xz(std::size_t a, std::size_t b, std::size_t i, std::size_t j, std::size_t k, std::size_t l);
};



} // namespace (nhfInt)