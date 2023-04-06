#include "nhfint/angmom_cart.hpp"
#include "nhfint/ericlass_cart.hpp"
#include "nhfint/int2e_repl_cart.hpp"
#include "nhfint/basis.hpp"
#include "nhfint/eritensor.hpp"
#include "nhfint/utility.hpp"
#include <cstddef>

namespace nhfInt {

EriTensorD0 int2e_repl_d0_cart(const SegmentBasisSet &bsSet) {
    EriTensorD0 eriD0(bsSet.cart_basis_size());

    for (std::size_t i = 0; i < bsSet.shell_size(); ++i) {
    for (std::size_t j = i; j < bsSet.shell_size(); ++j) {
        const std::size_t ij = idx2(i, j);
        const Shell &a = bsSet(i);
        const Shell &b = bsSet(j);
        
        for (std::size_t k = 0; k < bsSet.shell_size(); ++k) {
        for (std::size_t l = k; l < bsSet.shell_size(); ++l) {
            const std::size_t kl = idx2(k, l);
            const Shell &c = bsSet(k);
            const Shell &d = bsSet(l);

            if (kl > ij) {
                break;
            }

            // TODO
            // using some method calculate a EriClassD0Cart
            EriClassD0Cart ecD0(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);


            // TODO
            // using different methods to make calculation fast


            for (std::size_t ai = 0; ai < ecD0.nbsA; ++ai){
            for (std::size_t bj = 0; bj < ecD0.nbsB; ++bj){
                const std::size_t ab = idx2(a.bsBegCart + ai, 
                                            b.bsBegCart + bj);
                
                for (std::size_t ck = 0; ck < ecD0.nbsC; ++ck){
                for (std::size_t dl = 0; dl < ecD0.nbsD; ++dl){
                    const std::size_t cd = idx2(c.bsBegCart + ck,
                                                d.bsBegCart + dl);
                    
                    // position of this eri in a EriTensorD0
                    const std::size_t pos = idx2(ab, cd);
                    eriD0(pos) = ecD0(ai, bj, ck, dl);
                }}
            }}
        }}
    }}

    return eriD0;
}



EriTensorD1 int2e_repl_d1_cart(const SegmentBasisSet &bsSet) {
    EriTensorD1 etD1(bsSet.cart_basis_size());


    for (std::size_t i = 0; i < bsSet.shell_size(); ++i) {
    for (std::size_t j = i; j < bsSet.shell_size(); ++j) {
        const std::size_t ij = idx2(i, j);
        const Shell &a = bsSet(i);
        const Shell &b = bsSet(j);

        for (std::size_t k = 0; k < bsSet.shell_size(); ++k) {
        for (std::size_t l = k; l < bsSet.shell_size(); ++l) {
            const std::size_t kl = idx2(k, l);
            const Shell &c = bsSet(k);
            const Shell &d = bsSet(l);

            if (kl > ij) {
                break;
            }


            // TODO: some methods to calculate a EriClassD1Cart
            EriClassD1Cart ecD1(a.angMomTot, b.angMomTot, c.angMomTot, d.angMomTot);


            for (std::size_t ai = 0; ai < ecD1.nbsA; ++ai) {
            for (std::size_t bj = 0; bj < ecD1.nbsB; ++bj) {
            for (std::size_t ck = 0; ck < ecD1.nbsC; ++ck) {
            for (std::size_t dl = 0; dl < ecD1.nbsD; ++dl) {
                std::size_t idxI = a.bsBegCart + ai;
                std::size_t idxJ = b.bsBegCart + bj;
                std::size_t idxK = c.bsBegCart + ck;
                std::size_t idxL = d.bsBegCart + dl;


                etD1.d1_x(idxI, idxJ, idxK, idxL, 0) = ecD1.d1_x(ai, bj, ck, dl, 0);
                etD1.d1_x(idxJ, idxI, idxK, idxL, 0) = ecD1.d1_x(ai, bj, ck, dl, 1);
                etD1.d1_x(idxK, idxL, idxI, idxJ, 0) = ecD1.d1_x(ai, bj, ck, dl, 2);
                etD1.d1_x(idxL, idxK, idxI, idxJ, 0) = ecD1.d1_x(ai, bj, ck, dl, 3);

                etD1.d1_y(idxI, idxJ, idxK, idxL, 0) = ecD1.d1_y(ai, bj, ck, dl, 0);
                etD1.d1_y(idxJ, idxI, idxK, idxL, 0) = ecD1.d1_y(ai, bj, ck, dl, 1);
                etD1.d1_y(idxK, idxL, idxI, idxJ, 0) = ecD1.d1_y(ai, bj, ck, dl, 2);
                etD1.d1_y(idxL, idxK, idxI, idxJ, 0) = ecD1.d1_y(ai, bj, ck, dl, 3);

                etD1.d1_z(idxI, idxJ, idxK, idxL, 0) = ecD1.d1_z(ai, bj, ck, dl, 0);
                etD1.d1_z(idxJ, idxI, idxK, idxL, 0) = ecD1.d1_z(ai, bj, ck, dl, 1);
                etD1.d1_z(idxK, idxL, idxI, idxJ, 0) = ecD1.d1_z(ai, bj, ck, dl, 2);
                etD1.d1_z(idxL, idxK, idxI, idxJ, 0) = ecD1.d1_z(ai, bj, ck, dl, 3);
            }}}}

        }}
    }}

    return etD1;
}


}   // namespace (nhfInt)