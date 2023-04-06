#include "nhfint/basis.hpp"
#include "nhfint/calc_ericlass_dkr.hpp"
#include "nhfint/ericlass_cart.hpp"
#include <cstddef>
#include <cmath>
#include <vector>

namespace nhfInt {
namespace dkr {

EriClassD1Cart calc_ericlass_d1(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
) {
    EriClassD1Cart ecD1Cart(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);

    // resesrve space for roots and weights of rys polynomial
    std::size_t angMomSum = a.angMomTot + b.angMomTot + c.angMomTot + d.angMomTot;
    std::size_t N = (angMomSum + 1) / 2 + 1;    // N > (L + 1) / 2
    std::vector<double> RT(N+1);
    std::vector<double> WT(N+1);

    std::vector<ITensorD0> itD0x(N+1, ITensorD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD0> itD0y(N+1, ITensorD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD0> itD0z(N+1, ITensorD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD1> itD1x(N+1, ITensorD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD1> itD1y(N+1, ITensorD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));
    std::vector<ITensorD1> itD1z(N+1, ITensorD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot));


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

                    // write to ecD1
                    // loop for positon of derivative

                    // TODO:
                    // changing the order of loops
                    for (std::size_t q = 0; q < 4; ++q) {
                        double primEriValDx = 0.0;
                        double primEriValDy = 0.0;
                        double primEriValDz = 0.0;
                        for (std::size_t ni = 0; ni <= N; ++ni) {
                            primEriValDx += itD0y[ni](iy,jy,ky,ly) * itD0z[ni](iz,jz,kz,lz) * itD1x[ni](ix,jx,kx,lx, q) * WT[ni];
                            primEriValDy += itD0z[ni](iz,jz,kz,lz) * itD0x[ni](ix,jx,kx,lx) * itD1y[ni](iy,jy,ky,ly, q) * WT[ni];
                            primEriValDz += itD0x[ni](ix,jx,kx,lx) * itD0y[ni](iy,jy,ky,ly) * itD1z[ni](iz,jz,kz,lz, q) * WT[ni];
                        }
                        ecD1Cart.d1_x(i, j, k, l, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValDx;
                        ecD1Cart.d1_y(i, j, k, l, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValDy;
                        ecD1Cart.d1_z(i, j, k, l, q) += gsPdAB.Cab * gsPdCD.Cab * primEriValDz;
                    }

                }}
            }}
        }}
    }}
    
    return ecD1Cart;
}




} // namespace (dkr)
} // namespace (nhfInt)