#pragma once

#include "nhfint/basis.hpp"
#include "nhfint/ericlass_cart.hpp"
#include "nhfint/controller.hpp"


/* Taketa-Huzinaga-Oohata method for two-electron 
 * repulsion integrals. Ref: doi: 10.1143/JPSJ.21.2313
 *
 * Notes:
 * 1. The original version of `THO` method is to integrate each 
 *    primitive. After a simple modification, calculating a whole 
 *    eri class at the same time, the calculation speed can be 
 *    greatly improved.
 *
 * 2. The `tho` integral algorithm is very slow and we do not use 
 *    this algorithm in real calculations. The point of this 
 *    algorithm is to check the correctness of the code of 
 *    other integration algorithms.
 *
 * 3. Thanks to `Ivo Filot` for his guidance some years ago, 
 *    some of the code here refers to his repository:
 *    https://github.com/ifilot/hfcxx
 *    https://github.com/ifilot/pyqint
 */

namespace nhfInt {
namespace tho {

/*
 * segmented basis set
 */
EriClassD0Cart calc_ericlass_d0(
    const Shell &a,
    const Shell &b,
    const Shell &c,
    const Shell &d,
    const GaussSet &gsSet,
    const BasisSet &bsSet,
    const ShellSet &shSet,
    const GaussPairDataSet &gsPdSet,
    const BasisPairDataSet &bsPdSet,
    const ShellPairDataSet &shPdSet,
    const EriController &ctrl
);

EriClassD1Cart calc_ericlass_d1(
    const Shell &a, 
    const Shell &b, 
    const Shell &c, 
    const Shell &d, 
    const GaussSet &gsSet, 
    const BasisSet &bsSet, 
    const ShellSet &shSet, 
    const GaussPairDataSet &gsPdSet, 
    const BasisPairDataSet &bsPdSet, 
    const ShellPairDataSet &shPdSet, 
    const EriController &ctrl
);

EriClassD2Cart calc_ericlass_d2(
    const Shell &a, 
    const Shell &b, 
    const Shell &c, 
    const Shell &d, 
    const GaussSet &gsSet, 
    const BasisSet &bsSet, 
    const ShellSet &shSet, 
    const GaussPairDataSet &gsPdSet, 
    const BasisPairDataSet &bsPdSet, 
    const ShellPairDataSet &shPdSet, 
    const EriController &ctrl
);


// The coefficient of x^j in the expansion of (x+a)^l (x+b)^m
// ref: doi: 10.1143/JPSJ.21.2313 (eq 2.4)
double binomial_prefactor(
    std::size_t j, std::size_t l, std::size_t m, 
    double a, double b
);


}   // namespace (tho)
}   // namespace (nhfInt)