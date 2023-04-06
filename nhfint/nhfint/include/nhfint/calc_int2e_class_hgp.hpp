#pragma once

#include "nhfint/basis.hpp"
#include "nhfint/controller.hpp"
#include "nhfint/ericlass_cart.hpp"
#include <cstddef>
#include <vector>


namespace nhfInt {
namespace hgp {

/*
 * segmented basis set
 */
EriClassD0Cart calc_ericlass_d0(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const BasisPairDataSet &bsPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
);

EriClassD1Cart calc_ericlass_d1(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const BasisPairDataSet &bsPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
);

EriClassD2Cart calc_ericlass_d2(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const BasisPairDataSet &bsPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
);


// ====================================================================


// // use HGP method to calculate a eri class
// // doi: 10.1063/1.455553

void hgp_hrr_bra(
    EriClassD0Cart &abxx, 
    const std::vector<EriClassD0Cart> &e0xx, 
    nhfMath::Vec3d AB
);

void hgp_hrr_ket(
    EriClassD0Cart &xxcd, 
    const std::vector<EriClassD0Cart> &xxf0, 
    nhfMath::Vec3d CD
);


// single hgp hrr
void hgp_hrr(
    std::vector<double> &hrrOut,
    const std::vector<double> &hrrInp,
    double x, double y, double z,
    std::size_t a, std::size_t b
);


// generate [e0|f0]
void hgp_vrr_primitive(
    std::vector<double> &e0f0,
    std::size_t e, std::size_t f,
    double zetaAB, double kAB, 
    double Px, double Py, double Pz, 
    double PAx, double PAy, double PAz,
    double zetaCD, double kCD, 
    double Qx, double Qy, double Qz, 
    double QCx, double QCy, double QCz
);



}   // namespace (hgp)
}   // namespace (nhfInt)