#include "nhfint/calc_ericlass_tho.hpp"
#include "nhfint/angmom_cart.hpp"
#include "nhfint/basis.hpp"
#include "nhfint/ericlass_cart.hpp"
#include "nhfint/utility.hpp"
#include "boysfun/boysfun.hpp"
#include "nhfmath/mathfun.hpp"
#include <cmath>
#include <cstddef>
#include <vector>


namespace nhfInt {
namespace tho {

// EriClassD0Cart calc_ericlass_d0(
//     const Shell &a, 
//     const Shell &b, 
//     const Shell &c, 
//     const Shell &d, 
//     const GaussSet &gsSet, 
//     const BasisSet &bsSet, 
//     const ShellSet &shSet, 
//     const GaussPairDataSet &gsPdSet, 
//     const BasisPairDataSet &bsPdSet, 
//     const ShellPairDataSet &shPdSet, 
//     const EriController &ctrl
// ) {
//     assert(a.angMomTot <= AngMomCartMaxL);
//     assert(b.angMomTot <= AngMomCartMaxL);
//     assert(c.angMomTot <= AngMomCartMaxL);
//     assert(d.angMomTot <= AngMomCartMaxL);

//     EriClassD0Cart abcd(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);

//     // reserve space for boysfun
//     std::size_t angMomSum = a.angMomTot + b.angMomTot + c.angMomTot + d.angMomTot;
//     std::vector<double> Fm(angMomSum + 1);

//     for (std::size_t ai = 0; ai < a.constraction_size(); ++ai) {
//     for (std::size_t bj = 0; bj < b.constraction_size(); ++bj) {
//         const GaussPairData &gsPdAB = gsPdSet(a.gsBegCart + ai * a.cart_basis_size(),
//                                               b.gsBegCart + bj * b.cart_basis_size());

//         for (std::size_t ck = 0; ck < c.constraction_size(); ++ck) {
//         for (std::size_t dl = 0; dl < d.constraction_size(); ++dl) {
//             const GaussPairData &gsPdCD = gsPdSet(c.gsBegCart + ck * c.cart_basis_size(),
//                                                   d.gsBegCart + dl * d.cart_basis_size());

//             const double x = gsPdAB.zeta * gsPdCD.zeta
//                              / (gsPdAB.zeta + gsPdCD.zeta)
//                              * (gsPdAB.Pab - gsPdCD.Pab).len2();

//             // evaluate boys function
//             for (std::size_t n = 0; n <= angMomSum; ++n) {
//                 Fm[n] = boysFun::boysfun(n, x);
//             }
            
            
//             // calculate the contribution to each ERI
//             for (std::size_t i = 0; i < abcd.nbsA; ++i) {
//             for (std::size_t j = 0; j < abcd.nbsB; ++j) {
//                 const GaussPairData &gsPdij = gsPdSet(a.gsBegCart + ai * a.cart_basis_size() + i,
//                                                       b.gsBegCart + bj * b.cart_basis_size() + j);

//                 for (std::size_t k = 0; k < abcd.nbsC; ++k) {
//                 for (std::size_t l = 0; l < abcd.nbsD; ++l) {
//                     const GaussPairData &gsPdkl = gsPdSet(c.gsBegCart + ck * c.cart_basis_size() + k,
//                                                           d.gsBegCart + dl * d.cart_basis_size() + l);

//                     // calculate primitive ERI value
//                     double eriVal = 0.0;

//                     // TODO


//                     abcd(i, j, k, l) += gsPdij.Cab * gsPdkl.Cab * eriVal;
//                 }}
//             }}
//         }}
//     }}

//     return abcd;
// }


double binomial_prefactor(
    std::size_t j, std::size_t l, std::size_t m, 
    double a, double b
) {
    double ret = 0.0;
    for (std::size_t p = 0; p <= l && p <= j; ++p) {
        std::size_t q = j-p;
        if (q <= m) {
            ret += nhfMath::combination(l, p) *
                   nhfMath::combination(m, q) *
                   std::pow(a, l-p) * std::pow(b, m-q);
        }
    }
    return ret;
}


}   // namespace (tho)
}   // namespace (nhfInt)