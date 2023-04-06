#pragma once

#include "nhfmath/vec3d.hpp"
#include "nhfint/basisfile.hpp"
#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <string>

// TODO:
// some descriptions of these classes

/* 
 * 
 * 
 * 
 * 
 * 
 */ 

namespace nhfInt {

class Gauss;
class Basis;
class Shell;
class GaussSet;
class BasisSet;
class ShellSet;
class GaussPairData;
class BasisPairData;
class ShellPairData;
class GaussPairDataSet;
class BasisPairDataSet;
class ShellPairDataSet;
class SegmentBasisSet;
class GeneralBasisSet;


// primitive gaussian function (PGF)
class  Gauss {
public:
    double          alpha;
    double          coeff;
    AngMomCart      angMom;
    nhfMath::Vec3d  centre;

    Gauss();
    Gauss(double alpha, double coeff, 
          AngMomCart angMom, nhfMath::Vec3d centre);

    double norm() const;        // normalization constant
    
    // norm_alp * norm_lmn == norm
    // This processing is done to facilitate 
    // the implementation of the hgp algorithm.
    double norm_alp() const;    // alpha-part of the norm
    double norm_lmn() const;    // angMom-part of the norm
};


// constracted gaussian function (CGF)
class Basis {
public:
    std::size_t         atomIdx;
    AngMomCart          angMom;
    nhfMath::Vec3d      centre;
    std::vector<Gauss>  gsList;

    Basis();
    // Basis(const std::vector<Gauss> &gsList);

    double norm_lmn() const;    // angMom-part of the norm

    std::size_t gauss_size() const;

    // access operator
    const Gauss& operator()(std::size_t i) const;
};


// basis functions in same shell
class Shell {
public:
    std::size_t     atomIdx;            // atom index
    std::size_t     angMomTot;          // total angular momentum (L)
    std::size_t     gsBegCart;          // gauss func begin id (Cartesian)
    std::size_t     gsBegSphl;          // gauss func begin id (Spherical)
    std::size_t     bsBegCart;          // basis func begin id (Cartesian)
    std::size_t     bsBegSphl;          // basis func begin id (Spherical)
    std::size_t     cartBasisSize;      // (L+1)*(L+2) / 2
    std::size_t     sphlBasisSize;      // 2*L + 1
    std::size_t     cartGaussSize;      // alpha.size() * cartBasisSize
    std::size_t     sphlGaussSize;      // alpha.size() * sphlBasisSize
    
    nhfMath::Vec3d          centre;     // gaussian centre
    std::vector<double>     alpha;      // gaussian exponent
    std::vector<double>     coeff;      // Combination coefficient

    // TODO:
    // constructors

    std::size_t cart_basis_size() const;
    std::size_t sphl_basis_size() const;
    std::size_t cart_gauss_size() const;
    std::size_t sphl_gauss_size() const;
    std::size_t constraction_size() const;  // alpha.size()
};


// just a set of Gauss
class GaussSet {
public:
    std::vector<Gauss> gsList;

    GaussSet();
    GaussSet(const std::vector<Gauss> &gsList);

    std::size_t gauss_size() const;

    // access operator
    const Gauss& operator()(std::size_t i) const;
};


// just a set of Basis
class BasisSet {
public:
    std::vector<Basis> bsList;

    BasisSet();
    BasisSet(const std::vector<Basis> &bsList);

    std::size_t basis_size() const;

    // access operator
    const Basis& operator()(std::size_t i) const;
};


// just a set of Shell
class ShellSet {
public:
    std::vector<Shell> shList;

    ShellSet();
    ShellSet(const std::vector<Shell> &shList);

    std::size_t shell_size() const;

    // access operator
    const Shell& operator()(std::size_t i) const;
};


// pair data of two Gauss
// GaussPairData: record shell-pair data, see HGP-1988 eq 8, 9 and 15.
// Each pair of pimitive gaussians will have a PairData.
// Note:
// 1. Order-independent pair data will be recorded, that 
//    says PairData(a,b) and PairData(b,a) are the same.
// 2. Only gaussian centre and gaussian exponent of primitive
//    gausssian basis functions are needed to generate PairData.
class GaussPairData {
public:
    double              zeta;   // alpha1 + alpha2 (HGP eq 8)
    double              Kab;    // Kab (HGP eq 15)
    double              Aab;    // coeff1 * coeff2 * norm_alp_1 * norm_alp_2
    double              Bab;    // norm_lmn_1 * norm_lmn_2
    double              Cab;    // coeff1 * coeff2 * norm1 * norm2
    double              Qab;    // Cab^0.5 * [ab|ab]^0.5  used for screening
    nhfMath::Vec3d      Pab;    // (alpha1 * centre1 + alpha2 * centre2) / zeta

    GaussPairData();
    GaussPairData(const Gauss &gs1, const Gauss &gs2);
};

class GaussPairDataSet {
public:
    std::vector<GaussPairData> gsPdList;

    GaussPairDataSet(const GaussSet &gsSet);

    const GaussPairData& operator()(std::size_t ij) const;
    const GaussPairData& operator()(std::size_t i, std::size_t j) const;
};


// pair data of two Basis
class BasisPairData {
public:
    double      Qab;            // (ab|ab)^0.5 used for screening
    double      Bab;            // norm_lmn_1 * norm_lmn_2

    BasisPairData();
    BasisPairData(const Basis &bs1, const Basis &bs2);
};

class BasisPairDataSet {
public:
    std::vector<BasisPairData> bsPdList;

    BasisPairDataSet(const BasisSet &bsSet);

    const BasisPairData& operator()(std::size_t ij) const;
    const BasisPairData& operator()(std::size_t i, std::size_t j) const;
};


// pair data o two Shell
class ShellPairData {
public:
    double Qmx;     // max Qab

    ShellPairData();
    ShellPairData(const Shell &sh1, const Shell &sh2);
};

class ShellPairDataSet {
public:
    std::vector<ShellPairData> shPdList;

    ShellPairDataSet();
    ShellPairDataSet(const ShellSet &shSet);

    const ShellPairData& operator()(std::size_t ij) const;
    const ShellPairData& operator()(std::size_t i, std::size_t j) const;
};


// segmented basis set
class SegmentBasisSet {
public:
    std::vector<Shell>  shellList;

    /*        counters         */
    std::size_t     cartBasisSize;      // sum of Shell: (L+1)*(L+2) / 2
    std::size_t     sphlBasisSize;      // sum of Shell: 2*L + 1
    std::size_t     cartGaussSize;      // sum of Shell: alpha.size() * cartBasisSize
    std::size_t     sphlGaussSize;      // sum of Shell: alpha.size() * sphlBasisSize


    // TODO:
    // constructors

    // SegmentBasisSet();
    // SegmentBasisSet(std::size_t nBasis);

    // // Create a basis set for one atom.
    // SegmentBasisSet(
    //     std::string     basisFile,  // basis file
    //     std::string     atomType,   // atom type
    //     nhfMath::Vec3d  centre      // atom position
    // );

    // // Create a basis set for a molecule with several atoms
    // SegmentBasisSet(
    //     const std::vector<std::string>     &basisFileList,
    //     const std::vector<std::string>     &atomTypeList,
    //     const std::vector<nhfMath::Vec3d>  &centreList
    // );
    

    // access operator
    const Shell& operator()(std::size_t i) const;

    std::size_t cart_basis_size() const;
    std::size_t sphl_basis_size() const;
    std::size_t cart_gauss_size() const;
    std::size_t sphl_gauss_size() const;
    std::size_t shell_size() const;         // shellList.size()
};



// TODO

// One BasisBlock store informatiom of several shells.
// These shells have the same center, angular momentum
// and Gaussian exponent.
class Block {
public:

    /*
     * Common information in a block
     */
    std::size_t             gsBeg;      // gauss func begin id
    std::size_t             angMom;     // angular momentum
    nhfMath::Vec3d          centre;     // gaussian centre
    std::vector<double>     alpha;      // gaussian exponent


    // TODO:
    // merge the information in a struct
    std::vector<std::size_t>            bsBegCart;
    std::vector<std::size_t>            bsBegSphl;
    std::vector<std::vector<double>>    coeff;

    std::size_t cart_basis_size() const;    
    std::size_t sphl_basis_size() const;
    std::size_t gauss_size() const;
    std::size_t shell_size() const;
};


class GeneralBasisSet {
public:
    std::vector<Block>  blockList;

    GeneralBasisSet();
    GeneralBasisSet(std::size_t nBlock);

    // Create a basis set for one atom.
    GeneralBasisSet(
        std::string     basisFile,  // basis file
        std::string     atomType,   // atom type
        nhfMath::Vec3d  centre      // atom position
    );

    // Create a basis set for a molecule with several atoms
    GeneralBasisSet(
        const std::vector<std::string>     &basisFileList,
        const std::vector<std::string>     &atomTypeList,
        const std::vector<nhfMath::Vec3d>  &centreList
    );


    // access operator
    Block&  operator[](std::size_t i) const;

    std::size_t cart_basis_size() const;
    std::size_t sphl_basis_size() const;
    std::size_t cart_gauss_size() const;
    std::size_t sphl_gauss_size() const;
    std::size_t shell_size() const;
};


}   // namespace (nhfInt)
