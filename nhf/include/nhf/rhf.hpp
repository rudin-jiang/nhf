#pragma once


#include "nhf/molecule.hpp"
#include "nhfint/eritensor.hpp"
#include "nhfmath/matrix.hpp"

namespace nhf {

class RHF {
public:
    Molecule    mol;

    nhfMath::Matrix     ovlp;
    nhfMath::Matrix     kine;
    nhfMath::Matrix     nucl;
    nhfInt::EriTensorD0 repl;

    nhfMath::Matrix     C;  // 
    nhfMath::Matrix     P;  // density matrix

    bool check_convergence() const;
    void run(double dE);
    void opt(double dE);

    nhfMath::Matrix grad() const;   //
    nhfMath::Matrix hess() const;   //
};

} // namespace (nhf)


