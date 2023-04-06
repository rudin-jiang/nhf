#pragma once

#include <vector>
#include <cstddef>


namespace nhfInt {

class Int2eTensorD0 {
public:
    std::size_t nBasis;
    std::size_t pageSize;

    // s8 symmetry
    std::vector<double> e2D0; 

    Int2eTensorD0();
    Int2eTensorD0(std::size_t nBasis);

    double  d0(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d0(std::size_t i, std::size_t j, std::size_t k, std::size_t l);

    // TODO:
    // here may be some two-index or one-index codes.
};


class Int2eTensorD1 {
public:
    std::size_t nBasis;
    std::size_t pageSize;   // idx2(nBasis-1, nBasis-1) + 1
    std::size_t cubeSize;   // nBasis * pageSize
    std::size_t pileSize;   // nBasis * cubeSize

    std::vector<double> e2D1x;     // size = pileSize
    std::vector<double> e2D1y;     // size = pileSize
    std::vector<double> e2D1z;     // size = pileSize

    Int2eTensorD1();
    Int2eTensorD1(std::size_t nBasis);

    // derivative of first basis function
    double  d1x0(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d1x0(std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d1y0(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d1y0(std::size_t i, std::size_t j, std::size_t k, std::size_t l);
    double  d1z0(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;
    double& d1z0(std::size_t i, std::size_t j, std::size_t k, std::size_t l);

    // specify the basis function to be derived
    double  d1x(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1x(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1y(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1y(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l);
    double  d1z(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l) const;
    double& d1z(std::size_t a, std::size_t i, std::size_t j, std::size_t k ,std::size_t l);

private:
    std::size_t _d1_index(
        std::size_t a,
        std::size_t i, std::size_t j,
        std::size_t k, std::size_t l
    ) const;
};


class Int2eTensorD2 {
public:
    std::size_t nBasis;

    /*
     *  Data ordering:
     *
     */
    std::vector<double> e2D2xx;
    std::vector<double> e2D2yy;
    std::vector<double> e2D2zz;
    std::vector<double> e2D2xy;
    std::vector<double> e2D2yz;
    std::vector<double> e2D2zx;

    Int2eTensorD2();
    Int2eTensorD2(std::size_t nBasis);

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

private:
    std::size_t _d2_index(
        std::size_t a, std::size_t b, 
        std::size_t i, std::size_t j, 
        std::size_t k, std::size_t l
    ) const;
};

} // namespace (nhfInt)