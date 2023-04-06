#include "nhfint/basis.hpp"
#include "nhfint/calc_ericlass_hgp.hpp"
#include "nhfint/ericlass_cart.hpp"
#include "boysfun/boysfun.hpp"
#include "nhfmath/vec3d.hpp"
#include <cassert>
#include <cstddef>
#include <vector>

// static const double 

namespace nhfInt {
namespace hgp {

EriClassD0Cart calc_ericlass_d0(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const BasisPairDataSet &bsPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
) {
    // changed in call function
    // TODO:
    // make the function more general
    assert(a.angMomTot + b.angMomTot >= c.angMomTot + d.angMomTot);
    assert(a.angMomTot >= b.angMomTot);
    assert(c.angMomTot >= d.angMomTot);

    assert(a.angMomTot <= AngMomCartMaxL);
    assert(b.angMomTot <= AngMomCartMaxL);
    assert(c.angMomTot <= AngMomCartMaxL);
    assert(d.angMomTot <= AngMomCartMaxL);

    EriClassD0Cart abcd(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);
    std::size_t angMomSum = a.angMomTot + b.angMomTot + c.angMomTot + d.angMomTot;

    // reserve space for [ss|ss]^(m)
    std::vector<double> ssssm(angMomSum + 1);

    std::size_t eMin = a.angMomTot;
    std::size_t fMin = c.angMomTot;
    std::size_t eMax = a.angMomTot + b.angMomTot;
    std::size_t fMax = c.angMomTot + d.angMomTot; 

    // [e0|f0]  eMax: a.angMomTot + b.angMomTot
    //          fMax: c.angMomTot + d.angMomTot
    std::size_t nbsE = (eMax + 1) * (eMax + 2) / 2;
    std::size_t nbsF = (fMax + 1) * (fMax + 2) / 2;
    std::vector<double> primInt(nbsE * nbsF);

    // (e0|f0)  e: a.angMomTot, ..., a.angMomTot + b.angMomTot
    //          f: c.angMomTot, ..., c.angMomTot + d.angMomTot
    std::vector<std::vector<EriClassD0Cart>> 
    e0f0(b.angMomTot + 1, std::vector<EriClassD0Cart>(d.angMomTot + 1));

    for (std::size_t e = eMin, i = 0; e <= eMax; ++e, ++i) {
    for (std::size_t f = fMin, j = 0; f <= fMax; ++f, ++j) {
        e0f0[i][j] = EriClassD0Cart(e, 0, f, 0);
    }}

    // (e0|cd)  e: a.angMomTot, ..., a.angMomTot + b.angMomTot
    std::vector<EriClassD0Cart> e0cd(b.angMomTot + 1);
    for (std::size_t e = eMin, i = 0; e <= eMax; ++e, ++i) {
        e0cd[i] = EriClassD0Cart(e, 0, c.angMomTot, d.angMomTot);
    }

    // generate all constracted e0f0 needed
    for (std::size_t ai = 0; ai < a.constraction_size(); ++ai) {
    for (std::size_t bj = 0; bj < b.constraction_size(); ++bj) {
        const GaussPairData &gsPdAB = gsPdSet(a.gsBegCart + ai * a.cart_basis_size(),
                                              b.gsBegCart + bj * b.cart_basis_size());

        for (std::size_t ck = 0; ck < c.constraction_size(); ++ck) {
        for (std::size_t dl = 0; dl < d.constraction_size(); ++dl) {
            const GaussPairData &gsPdCD = gsPdSet(c.gsBegCart + ck * c.cart_basis_size(),
                                                  d.gsBegCart + dl * d.cart_basis_size());

            // coeffA * norm_alp_A * coeffB * norm_alp_B *
            // coeffC * norm_alp_C * coeffD * norm_alp_D
            const double coeffABCD = gsPdAB.Aab * gsPdCD.Aab;

            const double x = gsPdAB.zeta * gsPdCD.zeta
                             / (gsPdAB.zeta + gsPdCD.zeta)
                             * (gsPdAB.Pab - gsPdCD.Pab).len2();

            // evaluate boys function
            for (std::size_t n = 0; n <= angMomSum; ++n) {

                // TODO:
                // ssssm[n] = boysFun::boysfun(n, x);
            }

            for (std::size_t e = eMin, i = 0; e <= eMax; ++e, ++i) {
            for (std::size_t f = fMin, j = 0; f <= fMax; ++f, ++j) {
                
                // TODO:
                // calculate primitive integral e0f0 using hgp vrr
                // hgp_vrr(primInt, ssssm, ..., e, f);



                for (std::size_t ief = 0; ief < e0f0[e][f].size(); ++ief) {
                    e0f0[i][j](ief) += coeffABCD * primInt[ief];
                }
            }}
        }}
    }}


    // hgp_hrr_ket    
    for (std::size_t i = 0; i < e0cd.size(); ++i) {
        // TODO
        // hgp_hrr_ket(e0cd[i], e0f0[i], c.centre - d.centre);

        // TODO:
        // if (d == 0)

    }


    // TODO
    // hgp_hrr_bra
    // hgp_hrr_bra(abcd, e0cd, a.centre - b.centre);

    // TODO:
    // if (b == 0)


    // angular momentum related part
    for (std::size_t i = 0; i < abcd.nbsA; ++i) {
    for (std::size_t j = 0; j < abcd.nbsB; ++j) {
        const BasisPairData &bsPdij = bsPdSet(a.bsBegCart + i,
                                              b.bsBegCart + j);

        for (std::size_t k = 0; k < abcd.nbsC; ++k) {
        for (std::size_t l = 0; l < abcd.nbsD; ++l) {
            const BasisPairData &bsPdkl = bsPdSet(c.bsBegCart + k,
                                                  d.bsBegCart + l);

            // TODO:
            // change to: abcd(ij, kl) *= bsPdij.Bab * bsPdkl.Bab; ?
            abcd(i,j,k,l) *= bsPdij.Bab * bsPdkl.Bab;
        }}
    }}

    return abcd;
}




void hgp_hrr_bra(EriClassD0Cart &abxx, const std::vector<EriClassD0Cart> &e0xx, nhfMath::Vec3d AB) {

    std::size_t na = abxx.angA;
    std::size_t nb = abxx.angB;

    // TODO:
    // make the function more general
    assert(na >= nb);

    // check for angC and angD
    assert(abxx.angC == e0xx.front().angC);
    assert(abxx.angD == e0xx.front().angD);

    // check for e
    assert(e0xx.front().angA == na);
    assert(e0xx.back().angA == na + nb);
    assert(e0xx.front().angB == 0);
    assert(e0xx.back().angB == 0);
    
    // reserve space for inp
    // \sum_{i = na}^{na+nb}  (i+1)*(i+2) / 2
    // (1+nb) * (6+9*na+3*na*na+5*nb+3*na*nb+nb*nb) / 6;
    std::size_t hrrInpSize = (1+nb) * (6+9*na+3*na*na+5*nb+3*na*nb+nb*nb) / 6;
    std::vector<double> hrrInp(hrrInpSize);

    // reserve space for out
    std::size_t hrrOutSize = abxx.nbsA * abxx.nbsB;
    std::vector<double> hrrOut(hrrOutSize);

    std::size_t idx = 0;
    for (std::size_t k = 0; k < abxx.nbsC; ++k) {
    for (std::size_t l = 0; l < abxx.nbsD; ++l) {
        
        // TODO:
        // maybe using an index of kl ?

        // copy input
        idx = 0;
        for (std::size_t p = 0; p < e0xx.size(); ++p) {
            for (std::size_t i = 0; i < e0xx[p].nbsA; ++i) {
                hrrInp[idx++] = e0xx[p](i,0,k,l);
            }
        }

        // TODO:
        // hgp_hrr(hrrOut, hrrInp, AB.x, AB.y, AB.z, na, nb);

        // copy output
        idx = 0;
        for (std::size_t i = 0; i < abxx.nbsA; ++i) {
        for (std::size_t j = 0; j < abxx.nbsB; ++j) {
            abxx(i,j,k,l) = hrrOut[idx++];
        }}
    }}
}


void hgp_hrr_ket(EriClassD0Cart &xxcd, const std::vector<EriClassD0Cart> &xxf0, nhfMath::Vec3d CD) {
    std::size_t nc = xxcd.angC;
    std::size_t nd = xxcd.angD;

    // TODO:
    // make the function more general
    assert(nc >= nd);

    // check for angC and angD
    assert(xxcd.angA == xxf0.front().angA);
    assert(xxcd.angB == xxf0.front().angB);

    // check for e
    assert(xxf0.front().angC == nc);
    assert(xxf0.back().angC == nc + nd);
    assert(xxf0.front().angD == 0);
    assert(xxf0.back().angD == 0);

    // reserve space for inp
    // \sum_{i = nc}^{nc+nd}  (i+1)*(i+2) / 2
    // (1+nd) * (6+9*nc+3*nc*nc+5*nd+3*nc*nd+nd*nd) / 6;
    std::size_t hrrInpSize = (1+nd) * (6+9*nc+3*nc*nc+5*nd+3*nc*nd+nd*nd) / 6;
    std::vector<double> hrrInp(hrrInpSize);

    // reserve space for out
    std::size_t hrrOutSize = xxcd.nbsC * xxcd.nbsD;
    std::vector<double> hrrOut(hrrOutSize);

    std::size_t idx = 0;
    for (std::size_t i = 0; i < xxcd.nbsA; ++i) {
    for (std::size_t j = 0; j < xxcd.nbsB; ++j) {
        
        // TODO:
        // maybe using an index of ij ?

        // copy input
        idx = 0;
        for (std::size_t p = 0; p < xxf0.size(); ++p) {
            for (std::size_t k = 0; k < xxf0[p].nbsC; ++k) {
                hrrInp[idx++] = xxf0[p](i,j,k,0);
            }
        }

        // TODO:
        // hgp_hrr(hrrOut, hrrInp, CD.x, CD.y, CD.z, nc, nd);

        // copy output
        idx = 0;
        for (std::size_t k = 0; k < xxcd.nbsC; ++k) {
        for (std::size_t l = 0; l < xxcd.nbsD; ++l) {
            xxcd(i,j,k,l) = hrrOut[idx++];
        }}
    }}
}


void hgp_hrr_0_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_1_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_1_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_2_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_2_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_2_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_3_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_3_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_3_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_3_3_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_4_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_4_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_4_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_4_3_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_4_4_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_5_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_5_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_5_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_5_3_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_5_4_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_5_5_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_3_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_4_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_5_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_6_6_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_3_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_4_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_5_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_6_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_7_7_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_0_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_1_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_2_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_3_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_4_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_5_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_6_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_7_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);
void hgp_hrr_8_8_(std::vector<double> &hrrOut, const std::vector<double> &hrrInp, double x, double y, double z);


// single hgp hrr
void hgp_hrr(
    std::vector<double> &hrrOut,
    const std::vector<double> &hrrInp,
    double x, double y, double z,
    std::size_t a, std::size_t b
) {
    // TODO:
    // make the function more general
    assert(a >= b);

    assert(a <= 8);
    assert(b <= 8);
    assert(b != 0);

    // if (a == 0 && b == 0)   hgp_hrr_0_0_(hrrOut, hrrInp, x, y, z);
    // if (a == 1 && b == 0)   hgp_hrr_1_0_(hrrOut, hrrInp, x, y, z);
    if (a == 1 && b == 1)   hgp_hrr_1_1_(hrrOut, hrrInp, x, y, z);
    // if (a == 2 && b == 0)   hgp_hrr_2_0_(hrrOut, hrrInp, x, y, z);
    if (a == 2 && b == 1)   hgp_hrr_2_1_(hrrOut, hrrInp, x, y, z);
    if (a == 2 && b == 2)   hgp_hrr_2_2_(hrrOut, hrrInp, x, y, z);
    // if (a == 3 && b == 0)   hgp_hrr_3_0_(hrrOut, hrrInp, x, y, z);
    if (a == 3 && b == 1)   hgp_hrr_3_1_(hrrOut, hrrInp, x, y, z);
    if (a == 3 && b == 2)   hgp_hrr_3_2_(hrrOut, hrrInp, x, y, z);
    if (a == 3 && b == 3)   hgp_hrr_3_3_(hrrOut, hrrInp, x, y, z);
    // if (a == 4 && b == 0)   hgp_hrr_4_0_(hrrOut, hrrInp, x, y, z);
    if (a == 4 && b == 1)   hgp_hrr_4_1_(hrrOut, hrrInp, x, y, z);
    if (a == 4 && b == 2)   hgp_hrr_4_2_(hrrOut, hrrInp, x, y, z);
    if (a == 4 && b == 3)   hgp_hrr_4_3_(hrrOut, hrrInp, x, y, z);
    if (a == 4 && b == 4)   hgp_hrr_4_4_(hrrOut, hrrInp, x, y, z);
    // if (a == 5 && b == 0)   hgp_hrr_5_0_(hrrOut, hrrInp, x, y, z);
    if (a == 5 && b == 1)   hgp_hrr_5_1_(hrrOut, hrrInp, x, y, z);
    if (a == 5 && b == 2)   hgp_hrr_5_2_(hrrOut, hrrInp, x, y, z);
    if (a == 5 && b == 3)   hgp_hrr_5_3_(hrrOut, hrrInp, x, y, z);
    if (a == 5 && b == 4)   hgp_hrr_5_4_(hrrOut, hrrInp, x, y, z);
    if (a == 5 && b == 5)   hgp_hrr_5_5_(hrrOut, hrrInp, x, y, z);
    // if (a == 6 && b == 0)   hgp_hrr_6_0_(hrrOut, hrrInp, x, y, z);
    if (a == 6 && b == 1)   hgp_hrr_6_1_(hrrOut, hrrInp, x, y, z);
    if (a == 6 && b == 2)   hgp_hrr_6_2_(hrrOut, hrrInp, x, y, z);
    if (a == 6 && b == 3)   hgp_hrr_6_3_(hrrOut, hrrInp, x, y, z);
    if (a == 6 && b == 4)   hgp_hrr_6_4_(hrrOut, hrrInp, x, y, z);
    if (a == 6 && b == 5)   hgp_hrr_6_5_(hrrOut, hrrInp, x, y, z);
    if (a == 6 && b == 6)   hgp_hrr_6_6_(hrrOut, hrrInp, x, y, z);
    // if (a == 7 && b == 0)   hgp_hrr_7_0_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 1)   hgp_hrr_7_1_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 2)   hgp_hrr_7_2_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 3)   hgp_hrr_7_3_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 4)   hgp_hrr_7_4_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 5)   hgp_hrr_7_5_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 6)   hgp_hrr_7_6_(hrrOut, hrrInp, x, y, z);
    if (a == 7 && b == 7)   hgp_hrr_7_7_(hrrOut, hrrInp, x, y, z);
    // if (a == 8 && b == 0)   hgp_hrr_8_0_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 1)   hgp_hrr_8_1_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 2)   hgp_hrr_8_2_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 3)   hgp_hrr_8_3_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 4)   hgp_hrr_8_4_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 5)   hgp_hrr_8_5_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 6)   hgp_hrr_8_6_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 7)   hgp_hrr_8_7_(hrrOut, hrrInp, x, y, z);
    if (a == 8 && b == 8)   hgp_hrr_8_8_(hrrOut, hrrInp, x, y, z);
}


} // namespace (hgp)
} // namespace (nhfInt)