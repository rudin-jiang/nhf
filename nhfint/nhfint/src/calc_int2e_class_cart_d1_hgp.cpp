#include "nhfint/calc_ericlass_hgp.hpp"


namespace nhfInt {
namespace hgp {

EriClassD1Cart calc_ericlass_d1(
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

    EriClassD1Cart abcd(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);
    std::size_t angMomSum = a.angMomTot + b.angMomTot + c.angMomTot + d.angMomTot;

    // reserve space for [ss|ss]^(m)
    std::vector<double> ssssm(angMomSum + 2);

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



} // namespace (hgp)
} // namespace (nhfInt)