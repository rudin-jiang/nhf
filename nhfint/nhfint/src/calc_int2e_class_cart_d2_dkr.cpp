#include "nhfint/basis.hpp"
#include "nhfint/calc_ericlass_dkr.hpp"
#include "nhfint/ericlass_cart.hpp"
#include <cstddef>
#include <cmath>
#include <vector>

namespace nhfInt {
namespace dkr {


EriClassD2Cart calc_ericlass_d2(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
) {
    EriClassD2Cart ecD2Cart(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);

    // resesrve space for roots and weights of rys polynomial
    std::size_t angMomSum = a.angMomTot + b.angMomTot + c.angMomTot + d.angMomTot;
    std::size_t N = (angMomSum + 1) / 2 + 1;    // N > (L + 1) / 2
    std::vector<double> RT(N + 1);
    std::vector<double> WT(N + 1);

    std::vector<ITensorD0> itD0x(N+1, ITensorD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD0> itD0y(N+1, ITensorD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD0> itD0z(N+1, ITensorD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD1> itD1x(N+1, ITensorD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD1> itD1y(N+1, ITensorD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD1> itD1z(N+1, ITensorD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD2> itD2x(N+1, ITensorD2(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD2> itD2y(N+1, ITensorD2(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD2> itD2z(N+1, ITensorD2(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));


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


                    // write to ecD22
                    // loop for positon of derivative

                    // TODO:
                    // changing the order of loops

                    // xx, yy, zz
                    for (std::size_t p = 0; p < 4; ++p) {
                    for (std::size_t q = p; q < 4; ++q) {
                        double primEriValD2xx = 0.0;
                        double primEriValD2yy = 0.0;
                        double primEriValD2zz = 0.0;
                        for (std::size_t ni = 0; ni <= N; ++ni) {
                            primEriValD2xx += itD0y[ni](iy,jy,ky,ly) * itD0z[ni](iz,jz,kz,lz) * itD2x[ni](ix,jx,kx,lx, p, q) * WT[ni];
                            primEriValD2yy += itD0z[ni](iz,jz,kz,lz) * itD0x[ni](ix,jx,kx,lx) * itD2y[ni](iy,jy,ky,ly, p, q) * WT[ni];
                            primEriValD2zz += itD0x[ni](ix,jx,kx,lx) * itD0y[ni](iy,jy,ky,ly) * itD2z[ni](iz,jz,kz,lz, p, q) * WT[ni];
                        }

                        ecD2Cart.d2_xx(i, j, k, l, p, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValD2xx;
                        ecD2Cart.d2_yy(i, j, k, l, p, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValD2yy;
                        ecD2Cart.d2_zz(i, j, k, l, p, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValD2zz;
                    }}

                    // xy, yz, zx
                    for (std::size_t p = 0; p < 4; ++p) {
                    for (std::size_t q = 0; q < 4; ++q) {
                        double primEriValD2xy = 0.0;
                        double primEriValD2yz = 0.0;
                        double primEriValD2zx = 0.0;
                        for (std::size_t ni = 0; ni <= N; ++ni) {
                            primEriValD2xy += itD2x[ni](ix,jx,kx,lx, p, q) * itD2y[ni](iy,jy,ky,ly, p, q) * itD0z[ni](iz,jz,kz,lz) * WT[ni];
                            primEriValD2yz += itD2y[ni](iy,jy,ky,ly, p, q) * itD2z[ni](iz,jz,kz,lz, p, q) * itD0x[ni](ix,jx,kx,lx) * WT[ni];
                            primEriValD2zx += itD2z[ni](iz,jz,kz,lz, p, q) * itD2x[ni](ix,jx,kx,lx, p, q) * itD0y[ni](iy,jy,ky,ly) * WT[ni];
                        }

                        ecD2Cart.d2_xy(i, j, k, l, p, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValD2xy;
                        ecD2Cart.d2_yz(i, j, k, l, p, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValD2yz;
                        ecD2Cart.d2_zx(i, j, k, l, p, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValD2zx;
                    }}

                }}
            }}
        }}
    }}
    
    return ecD2Cart;
}


} // namespace (dkr)
} // namespace (nhfInt)