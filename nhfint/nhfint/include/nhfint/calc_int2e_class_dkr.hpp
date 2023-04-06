#pragma once

#include "nhfint/basis.hpp"
#include "nhfint/controller.hpp"
#include "nhfint/ericlass_cart.hpp"
#include <cstddef>
#include <vector>

namespace nhfInt {
namespace dkr {

/*
 * segmented basis set
 */
EriClassD0Cart calc_ericlass_d0(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
);

EriClassD1Cart calc_ericlass_d1(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
);

EriClassD2Cart calc_ericlass_d2(
    const Shell &a, const Shell &b,
    const Shell &c, const Shell &d,
    const ShellPairDataSet &shPdSet,
    const GaussPairDataSet &gsPdSet,
    const EriController &ctrl
);

// TODO: testing
// TODO: better comment
// In a EriClass, all terms I needed are like 
// a 4D tensor.
// I(i,j,k,l)
// i range from 0 to angA,  j range from 0 to angB
// k range from 0 to angC,  l range from 0 to angD
// attention: i can equal angA, ...
// ref: 10.1002/jcc.540110809
class ITD0 {
public:
    std::size_t angA, angB, angC, angD;
    std::size_t cntA, cntB, cntC, cntD; // cntX = angX + 1
    std::size_t num0, num1, num2, num3; // fasten indexing
    std::vector<double> itD0;

    ITD0();
    ITD0(std::size_t angA, std::size_t angB, 
         std::size_t angC, std::size_t angD);

    // view ITensor as four-tuple
    // i range from 0 to cntA - 1
    // j range from 0 to cntB - 1
    // k range from 0 to cntC - 1
    // l range from 0 to cntD - 1
    double  operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l) const;
    double& operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l);

    // view ITensor as two-tuple
    // ij range from 0 to cntA * cntB - 1
    // kl range from 0 to cntC * cntD - 1
    double  operator()(std::size_t ij, std::size_t kl) const;
    double& operator()(std::size_t ij, std::size_t kl);

    // get data directly
    // ijkl range from 0 to cntA * cntB * cntC * cntD - 1
    double  operator()(std::size_t ijkl) const;
    double& operator()(std::size_t ijkl);
};



class ITD1 {
public:
    std::size_t angA, angB, angC, angD;
    std::size_t cntA, cntB, cntC, cntD; // cntX = angX + 1
    std::size_t num0, num1, num2, num3; // fasten indexing
    std::vector<double> itD1;

    ITD1();
    ITD1(std::size_t angA, std::size_t angB, 
         std::size_t angC, std::size_t angD);

    double  operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l, 
                                      std::size_t a) const;
    double& operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l,
                                      std::size_t a);
};


class ITD2 {
public:
    std::size_t angA, angB, angC, angD;
    std::size_t cntA, cntB, cntC, cntD; // 
    std::size_t num0, num1, num2, num3; // fasten indexing
    std::vector<double> itD2;

    ITD2();
    ITD2(std::size_t angA, std::size_t angB, 
         std::size_t angC, std::size_t angD);

    double  operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l,
                       std::size_t a, std::size_t b) const;
    double& operator()(std::size_t i, std::size_t j,
                       std::size_t k, std::size_t l,
                       std::size_t a, std::size_t b);
};


// transform ITensor(na+nb, 0, x, x) to ITenosr(na, nb, x, x)
// ref: 10.1002/jcc.540110809 (eq 19)
// Note: eq 19 is error
//  qn = (jx, n) * (xi - xj) ^ (jx - n)
ITD0 dkr_hrr_bra(const ITD0 &n0xx, double Xij, std::size_t na, std::size_t nb);

// transform ITensor(x, x, nc+nd, 0) to ITensor(x, x, nc, nd)
// ref: 10.1002/jcc.540110809 (like eq 19)
ITD0 dkr_hrr_ket(const ITD0 &xxm0, double Xkl, std::size_t nc, std::size_t nd);

// TODO: test
// transform ITensor(na+nb, 0, nc+nd, 0) to ITensor(na, nb, nc, nd)
ITD0 dkr_hrr(const ITD0 &n0m0, double Xji, double Xlk, 
             std::size_t na, std::size_t nb, std::size_t nc, std::size_t nd);

// generate n0m0
// for all n ranges from 0 to nMax,
// for all m ranges from 0 to mMax
// ref: 10.1002/jcc.540110809 (eq15 and eq16)
// Here I use a new name that is different from the original ref.
// I0 = I(0,0,0,0);  B0 = B00;  BB = B10;  BP = B10';  CC = C00;  CP = C00';
ITD0 dkr_vrr(double I0, double B0, double BB, double BP,
             double CC, double CP, std::size_t nMax, std::size_t mMax);

}   // namespace (dkr)
}   // namespace (nhfInt)