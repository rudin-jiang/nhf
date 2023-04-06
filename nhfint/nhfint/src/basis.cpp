#include "nhfint/basis.hpp"
#include "nhfint/angmom_cart.hpp"
#include "nhfmath/vec3d.hpp"
#include "nhfint/utility.hpp"
#include <cassert>
#include <vector>
#include <string>
#include <cstddef>
#include <cmath>
#include <fstream>

// // pow(2.0, 0.5) * pow(pi, 1.25), used in KAB and KCD
// // in function hgp_vrr_contracted
// static const double preCoeff = 5.9149671727956132;

namespace nhfInt {



// std::size_t Shell::cart_basis_size() const
// { return (angMom+1) * (angMom+2) / 2; }

// std::size_t Shell::sphl_basis_size() const
// { return 2 * angMom + 1; }

// std::size_t Shell::gauss_size() const
// { return alpha.size(); }


// Basis::Basis()
// : bsBeg(0), gsBeg(0), angMom(0) {}

// Basis::Basis(
//     std::size_t bsBeg, std::size_t gsBeg, std::size_t angMom,
//     const VecReal &alpha, const VecReal &coeff, const Vec3d &centre
// ) : bsBeg(bsBeg), gsBeg(gsBeg), angMom(angMom),
//     alpha(alpha), coeff(coeff), centre(centre) {}


// /*          BasisSet            */
// BasisSet::BasisSet() {}

// BasisSet::BasisSet(std::size_t nBasis)
// : bsData(std::vector<Basis>(nBasis)) {}

// // BasisSet::BasisSet(std::string basisFile, std::string atomType, const Vec3d &centre) {}


// Basis   BasisSet::operator[](std::size_t i) const
// { return bsData[i]; }

// Basis&  BasisSet::operator[](std::size_t i)
// { return bsData[i]; }

// std::size_t BasisSet::basis_size() const{
//     std::size_t nBs = 0;
//     for (const Basis &b : bsData)
//         nBs += b.basis_size();
//     return nBs;
// }

// std::size_t BasisSet::gauss_size() const {
//     std::size_t nGs = 0;
//     for (const Basis &b : bsData)
//         nGs += b.gauss_size();
//     return nGs;
// }

// std::size_t BasisSet::shell_size() const 
// { return bsData.size(); }


// /*          PairData            */
// PairData::PairData()
// : zeta(0.0), Kab(0.0), Cab(0.0) {}

// PairData::PairData(double zeta, double Kab, double Cab, const Vec3d &P)
// : zeta(zeta), Kab(Kab), Cab(Cab), P(P) {}


// /*          PairDataSet         */
// PairDataSet::PairDataSet() {}

// PairDataSet::PairDataSet(std::size_t nPairData)
// : pdData(std::vector<PairData>(nPairData)) {}

// PairDataSet::PairDataSet(const BasisSet &bss) {
//     std::size_t nGs = bss.gauss_size();
//     std::size_t nPr = idx2(nGs-1, nGs-1) + 1;
//     pdData = std::vector<PairData>(nPr);

//     for (std::size_t i = 0; i < bss.shell_size(); ++i) {
//         const Basis &a = bss[i];

//     for (std::size_t j = 0; j <= i; ++j) {
//         const Basis &b = bss[j];

//         for (std::size_t ia = 0; ia < a.gauss_size(); ++ia) {
//         for (std::size_t jb = 0; jb < b.gauss_size(); ++jb) {

//             double zeta = a.alpha[ia] + b.alpha[jb];
//             double invz = 1.0 / zeta;
//             double len2 = (a.centre - b.centre).len2();
//             double inp  = -a.alpha[ia] * b.alpha[jb] * invz * len2;
//             double Kab  = preCoeff * invz * std::exp(inp);
//             double Cab  = a.coeff[ia] * b.coeff[jb];
//             Vec3d  P    = (a.alpha[ia] * a.centre + b.alpha[jb] * b.centre) * invz;

//             // position of this pair in PairDataSet
//             std::size_t iPos = idx2(a.gsBeg + ia, b.gsBeg + jb);
//             pdData[iPos] = PairData(zeta, Kab, Cab, P);
//         }}

//     }}
// }

// PairData PairDataSet::operator[](std::size_t i) const
// { return pdData[i]; }

// PairData& PairDataSet::operator[](std::size_t i)
// { return pdData[i]; }



// Gauss::Gauss(): alpha(0.0), coeff(0.0) {}

// Gauss::Gauss(double alpha, double coeff, AngMomCart angMom, nhfMath::Vec3d centre)
// : alpha(alpha), coeff(coeff), angMom(angMom), centre(centre) {
//     assert(alpha > 0.0);
// }

// double Gauss::norm() const {
//     return gauss_norm_cart(alpha, angMom.lx, angMom.ly, angMom.lz);
// }


// GaussSet::GaussSet(const SegmentBasisSet &segBsSet) {
    
// }

// const Gauss& GaussSet::operator()(std::size_t i) const {
//     assert(i < gsList.size());
//     return gsList[i];
// }






}   // namespace (nhfInt)