#include "nhfint/calc_ericlass_dkr.hpp"
#include "nhfint/ericlass_cart.hpp"
#include "nhfint/angmom_cart.hpp"
#include "nhfint/basis.hpp"
#include "nhfmath/mathfun.hpp"
#include <cstddef>
#include <vector>
#include <cassert>

namespace nhfInt {
namespace dkr {

EriClassD0Cart calc_ericlass_d0(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
) {
    EriClassD0Cart ecD0Cart(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);

    // resesrve space for roots and weights of rys polynomial
    std::size_t angMomSum = a.angMomTot + b.angMomTot + c.angMomTot + d.angMomTot;
    std::size_t N = (angMomSum + 1) / 2 + 1;    // N > (L + 1) / 2
    std::vector<double> RT(N+1);
    std::vector<double> WT(N+1);

    std::vector<ITD0> itD0x(N+1, ITD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITD0> itD0y(N+1, ITD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITD0> itD0z(N+1, ITD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));

    for (std::size_t ai = 0; ai < a.constraction_size(); ++ai) {
    for (std::size_t bj = 0; bj < b.constraction_size(); ++bj) {
        const ShellPairData &shPdAB = shPdSet(a.gsBegCart + ai,
                                              b.gsBegCart + bj);

        for (std::size_t ck = 0; ck < c.constraction_size(); ++ck) {
        for (std::size_t dl = 0; dl < d.constraction_size(); ++dl) {
            const ShellPairData &shPdCD = shPdSet(c.gsBegCart + ck,
                                                  d.gsBegCart + dl);

            // TODO:
            // screening
            // double EST = shPdAB.Qmx * shPdCD.Qmx;
            // if (ctrl.eriScreenSwitch && EST < ctrl.eriScreenCutoff) {
            //     continue;
            // }

            // TODO:
            // calculate roots and weights of rys polynomial


            // TODO:
            // 1. calculate itD0 and itD1
            // 2. Reduce repeated allocation of memory



            for (std::size_t i = 0; i < a.cart_basis_size(); ++i) {
            for (std::size_t j = 0; j < b.cart_basis_size(); ++j) {
                const GaussPairData &gsPdAB = gsPdSet(a.gsBegCart + ai * a.cart_basis_size() + i,
                                                      b.gsBegCart + bj * b.cart_basis_size() + j);

                for (std::size_t k = 0; k < c.cart_basis_size(); ++k) {
                for (std::size_t l = 0; l < d.cart_basis_size(); ++l) {
                    const GaussPairData &gsPdCD = gsPdSet(c.gsBegCart + ck * c.cart_basis_size() + k,
                                                          d.gsBegCart + dl * d.cart_basis_size() + l);

                    // write to ecD1
                    // loop for positon of derivative

                    // TODO:

                    std::size_t ix = angmomcart_component_x(a.angMomTot, i);
                    std::size_t iy = angmomcart_component_y(a.angMomTot, i);
                    std::size_t iz = angmomcart_component_z(a.angMomTot, i);
                    std::size_t jx = angmomcart_component_x(b.angMomTot, j);
                    std::size_t jy = angmomcart_component_y(b.angMomTot, j);
                    std::size_t jz = angmomcart_component_z(b.angMomTot, j);
                    std::size_t kx = angmomcart_component_x(c.angMomTot, k);
                    std::size_t ky = angmomcart_component_y(c.angMomTot, k);
                    std::size_t kz = angmomcart_component_z(c.angMomTot, k);
                    std::size_t lx = angmomcart_component_x(d.angMomTot, l);
                    std::size_t ly = angmomcart_component_y(d.angMomTot, l);
                    std::size_t lz = angmomcart_component_z(d.angMomTot, l);


                    double primEriVal = 0.0;
                    for (std::size_t ni = 0; ni <= N; ++ni) {
                        primEriVal += itD0x[ni](ix,jx,kx,lx) * itD0y[ni](iy,jy,ky,ly) * itD0z[ni](iz,jz,kz,lz) * WT[ni];
                    }

                    ecD0Cart(i, j, k, l) += gsPdAB.Cab * gsPdCD.Cab * primEriVal;
                }}
            }}
        }}
    }}

    return ecD0Cart;
}


// TODO: testing
ITD0::ITD0(): 
angA(0), angB(0), angC(0), angD(0),
cntA(1), cntB(1), cntC(1), cntD(1),
num0(1), num1(1), num2(1), num3(1),
itD0(num0, 0.0) {}

ITD0::ITD0(std::size_t angA, std::size_t angB, 
                 std::size_t angC, std::size_t angD)
: angA(angA), angB(angB), angC(angC), angD(angD),
  cntA(angA+1), cntB(angB+1), cntC(angC+1), cntD(angD+1)
{   
    num3 = cntD;
    num2 = cntC * num3;
    num1 = cntB * num2;
    num0 = cntA * num1;
    itD0 = std::vector<double>(num0, 0.0);
}

double ITD0::operator()(std::size_t i, std::size_t j,
                           std::size_t k, std::size_t l) const {
    assert(i < cntA);
    assert(j < cntB);
    assert(k < cntC);
    assert(l < cntD);
    return itD0[i * num1 + j * num2 + k * num3 + l];
}

double& ITD0::operator()(std::size_t i, std::size_t j,
                            std::size_t k, std::size_t l) {
    assert(i < cntA);
    assert(j < cntB);
    assert(k < cntC);
    assert(l < cntD);
    return itD0[i * num1 + j * num2 + k * num3 + l];
}

double  ITD0::operator()(std::size_t ij, std::size_t kl) const {
    assert(ij < cntA * cntB);
    assert(kl < cntC * cntD);
    return itD0[ij * num2 + kl];
}

double& ITD0::operator()(std::size_t ij, std::size_t kl) {
    assert(ij < cntA * cntB);
    assert(kl < cntC * cntD);
    return itD0[ij * num2 + kl];
}

double  ITD0::operator()(std::size_t ijkl) const {
    assert(ijkl < cntA * cntB * cntC * cntD);
    return itD0[ijkl];
}

double& ITD0::operator()(std::size_t ijkl) {
    assert(ijkl < cntA * cntB * cntC * cntD);
    return itD0[ijkl];
}



ITD0 dkr_hrr_bra(const ITD0 &n0xx, double Xij, std::size_t na, std::size_t nb) {
    assert(n0xx.angA == na + nb);
    assert(n0xx.angB == 0);

    // TODO: change to: double xPow[10] ?
    // calculate pow(Xji, n), n from 0 to nb
    std::vector<double> xPow(nb+1, 1.0);
    for (std::size_t i = 1; i < xPow.size(); ++i) {
        xPow[i] = xPow[i-1] * Xij;
    }

    ITD0 abxx(na, nb, n0xx.angC, n0xx.angD);
    for (std::size_t jx = 0; jx < abxx.cntB; ++jx) {
        for (std::size_t n = 0; n <= jx; ++n) {
            double qn = nhfMath::combination(jx, n) * xPow[jx - n];
            for (std::size_t ix = 0; ix < abxx.cntA; ++ix) {
                std::size_t ij_n0xx = (ix + n) * n0xx.num1;
                std::size_t ij_abxx = ix * abxx.num1 + jx * abxx.num2;

                for (std::size_t kl = 0; kl < abxx.num2; ++kl) {
                    abxx(ij_abxx + kl) += qn * n0xx(ij_n0xx + kl);
                }
            }
        }
    }

    return abxx;
}


ITD0 dkr_hrr_ket(const ITD0 &xxm0, double Xkl, std::size_t nc, std::size_t nd) {
    assert(xxm0.angC == nc + nd);
    assert(xxm0.angD == 0);

    // TODO: change to: double xPow[10] ?
    // calculate pow(Xlk, n), n from 0 to nd
    std::vector<double> xPow(nd+1, 1.0);
    for (std::size_t i = 1; i < xPow.size(); ++i) {
        xPow[i] = xPow[i-1] * Xkl;
    }

    ITD0 xxcd(xxm0.angA, xxm0.angB, nc, nd);
    for (std::size_t lx = 0; lx < xxcd.cntD; ++lx) {
        for (std::size_t n = 0; n <= lx; ++n) {
            double qn = nhfMath::combination(lx, n) * xPow[lx - n];

            for (std::size_t kx = 0; kx < xxcd.cntC; ++kx) {
                std::size_t kl_xxm0 = (kx + n) * xxm0.num3;
                std::size_t kl_xxcd = kx * xxcd.num3 + lx;

                for (std::size_t ij = 0; ij < xxcd.cntA * xxcd.cntB; ++ij) {
                    xxcd(ij * xxcd.num2 + kl_xxcd) += qn * xxm0(ij * xxm0.num2 + kl_xxm0);
                }
            }
        }
    }

    return xxcd;
}


ITD0 dkr_hrr(const ITD0 &n0m0, double Xij, double Xkl, 
                std::size_t na, std::size_t nb, std::size_t nc, std::size_t nd) {
    assert(n0m0.angA == na + nb);
    assert(n0m0.angC == nc + nd);
    assert(n0m0.angB == 0);
    assert(n0m0.angD == 0);

    ITD0 abm0 = dkr_hrr_bra(n0m0, Xij, na, nb);
    ITD0 abcd = dkr_hrr_ket(abm0, Xkl, nc, nd);

    return abcd;
}


ITD0 dkr_vrr(double I0, double B0, double BB, double BP,
                double CC, double CP, std::size_t n, std::size_t m) {
    ITD0 n0m0(n, 0, m, 0);

    // initial value
    n0m0(0, 0) = I0;

    // for (std::size_t )
}




}   // namespace (dkr)
}   // namespace (nhfInt)