
#include "cg_base/ijk.hpp"
#include "cg_base/iterm.hpp"
#include <cstddef>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>

std::vector<ITerm> iterm_needed_by_one_ericlass(
    std::size_t na, std::size_t nb, std::size_t nc, std::size_t nd
) {
    std::vector<ijk> shellA = generate_shell(na);
    std::vector<ijk> shellB = generate_shell(nb);
    std::vector<ijk> shellC = generate_shell(nc);
    std::vector<ijk> shellD = generate_shell(nd);

    std::unordered_set<ITerm, iterm_hash> itermSet;
    for (const ijk &a : shellA) {
    for (const ijk &b : shellB) {
    for (const ijk &c : shellC) {
    for (const ijk &d : shellD) {
        itermSet.insert(ITerm(a.i, b.i, c.i, d.i));
    }}}}

    std::vector<ITerm> itermVec(itermSet.begin(), itermSet.end());
    std::sort(itermVec.begin(), itermVec.end());
    
    return itermVec;
}


int main() {
    const std::size_t angMomMax = 6;
    
    // guess that: in an ericlass (na, nb, nc, nd)
    // all ITerms (i, j, k, l) needed are Iterms that
    // i range from 0 to na,  j range from 0 to nb,
    // k range from 0 to nc,  l range from 0 to nd.
    for (std::size_t na = 0; na <= angMomMax; ++na){
    for (std::size_t nb = 0; nb <= angMomMax; ++nb){
    for (std::size_t nc = 0; nc <= angMomMax; ++nc){
    for (std::size_t nd = 0; nd <= angMomMax; ++nd){
        std::vector<ITerm> itermNeedBySimul = iterm_needed_by_one_ericlass(na, nb, nc, nd);
        std::vector<ITerm> itermNeedByGuess = generate_iterm(na, nb, nc, nd);

        if (itermNeedBySimul != itermNeedByGuess) {
            std::printf("Guess not established in EriClass (%zu, %zu, %zu, %zu)\n",
                                                                        na, nb, nc, nd);
        }
    }}}}

    std::cout << "Hello, world\n" <<std::endl;
    return 0;
}





