#include "cg_base/ijkpair.hpp"
#include "cg_base/ijkpairm.hpp"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <random>
#include <fstream>
#include <unordered_set>
#include <vector>
#include <queue>

std::vector<ijkPairM> hgpVrrPathMin;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<std::size_t> dis(0);

std::vector<ijkPairM> hgp_vrr_need(const ijkPairM &pm, std::size_t i) {
    assert(i < 6);
    assert(pm[i] != 0);

    std::size_t j = (i+3) % 6;

    std::vector<ijkPairM> ret;

    ijkPairM tmp1 = pm;
    tmp1[i] -= 1;
    ret.push_back(tmp1);
    tmp1.m += 1;
    ret.push_back(tmp1);

    if (pm[i] >= 2) {
        ijkPairM tmp2 = pm;
        tmp2[i] -= 2;
        ret.push_back(tmp2);
        tmp2.m += 1;
        ret.push_back(tmp2);
    }

    if (pm[j] != 0) {
        tmp1[j] -= 1;
        ret.push_back(tmp1);
    }

    return ret;
}


std::vector<ijkPairM> hgp_vrr_path_gen_random(std::size_t na, std::size_t nb) {
    assert(na <= 16);
    assert(nb <= 16);
    assert(na >= nb);

    std::unordered_set<ijkPairM, ijkPairM_hash> vrrPath;
    std::unordered_set<ijkPairM, ijkPairM_hash> waitForGen;

    std::vector<ijkPairM> vrrNeed = generate_ijkpair0(na, nb);

    // 2 zeros or 3 zeros
    for (const ijkPairM &pm : vrrNeed) {
        waitForGen.insert(pm);
    }

    while (! waitForGen.empty()) {
        ijkPairM tmp = *(waitForGen.begin());
        waitForGen.erase(tmp);

        // insert to path first
        vrrPath.insert(tmp);

        if (tmp.cnt_nonz() == 0) {
            continue;
        }

        std::vector<std::size_t> nonzIdx = tmp.nonzero_index();
        std::size_t dfsIdx = nonzIdx[dis(gen) % nonzIdx.size()];
        std::vector<ijkPairM> ab = hgp_vrr_need(tmp, dfsIdx);
        for (const ijkPairM &pm : ab) {
            waitForGen.insert(pm);
        }
    }

    std::vector<ijkPairM> pairmVec(vrrPath.begin(), vrrPath.end());
    return pairmVec;
}


int main(int argc, char **argv) {
    
    if (argc != 3) {
        std::cerr << "error argumnets." << std::endl;
        return -1;
    }
    
    std::size_t nTest = 1000;
    
    std::size_t na = std::stoi(argv[1]);
    std::size_t nb = std::stoi(argv[2]);

    if (na > 16 || nb > 16 || nb > na) {
        std::cerr << "error argumnets." << std::endl;
        return -1;
    }

    // for (std::size_t n = 0; n < nTest; ++n) {
    //     std::printf("na = %2zu    nb = %2zu    try case: %5zu / %5zu\n",
    //                                                     na, nb, n, nTest);
    //     std::vector<ijkPair> pathThisTest = hgp_hrr_path_gen_random(na, nb);

    //     if (n == 0 || pathThisTest.size() < hgpHrrPathMin.size()) {
    //         hgpHrrPathMin = pathThisTest;
    //     }
    // }

    hgpVrrPathMin = hgp_vrr_path_gen_random(na, nb);



    // std::sort(hgpHrrPathMin.begin(), hgpHrrPathMin.end(), [](const ijkPair &pa, const ijkPair &pb){
    //     if (pa.b.cnt_zero() != pb.b.cnt_zero()) return pa.b.cnt_zero() < pb.b.cnt_zero();
    //     return pa < pb;
    // });

    std::string filename = "hgp_vrr_path_" + std::to_string(na) + "_" + std::to_string(nb) + 
                            "_" + std::to_string(hgpVrrPathMin.size()) + ".txt";
    std::ofstream ofs(filename);
    for (const ijkPairM &pm : hgpVrrPathMin) {
        ofs << pm.to_string() << std::endl;
    }
    ofs.close();

    return 0;
}