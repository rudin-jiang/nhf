#include "cg_base/ijkpair.hpp"
#include "cg_base/ijk.hpp"
#include <cassert>
#include <cstddef>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <queue>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>

std::vector<ijkPair> hgpHrrPathMin;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<std::size_t> dis(0);

// pos range from 0 to 2
std::vector<ijkPair> hgp_hrr_need(const ijkPair &p, std::size_t pos) {
    assert(p.b[pos] != 0);

    ijkPair needA = p;
    ijkPair needB = p;

    needA.a[pos] += 1;
    needA.b[pos] -= 1;
    needB.b[pos] -= 1;

    return {needA, needB};
}

std::vector<ijkPair> hgp_hrr_path_gen_random(std::size_t na, std::size_t nb) {
    assert(na <= 8);
    assert(nb <= 8);
    assert(na >= nb);

    std::unordered_set<ijkPair, ijkPair_hash> hrrPath;
    std::unordered_set<ijkPair, ijkPair_hash> pairSet;

    std::vector<ijkPair> hrrNeed = generate_pair(na, nb);
    
    // 0 zero or 1 zero
    for (const ijkPair &p : hrrNeed) {
        if (p.b.cnt_nonz() <= 1) {
            pairSet.insert(p);
        }
    }

    while (! pairSet.empty()) {
        ijkPair tmp = *(pairSet.begin());
        pairSet.erase(tmp);

        // insert to path first
        hrrPath.insert(tmp);

        if (tmp.b.cnt_nonz() == 0) {
            continue;
        }

        assert(tmp.b.cnt_nonz() == 1);
        std::size_t dfsIdx = tmp.b.first_nonzero_index();
        std::vector<ijkPair> ab = hgp_hrr_need(tmp, dfsIdx);
        for (const ijkPair &p : ab) {
            pairSet.insert(p);
        }
    }


    // 2 zeros or 3 zeros
    for (const ijkPair &p : hrrNeed) {
        if (p.b.cnt_nonz() >= 2) {
            pairSet.insert(p);
        }
    }

    while (! pairSet.empty()) {
        ijkPair tmp = *(pairSet.begin());
        pairSet.erase(tmp);

        // insert to path first
        hrrPath.insert(tmp);

        if (tmp.b.cnt_nonz() == 0) {
            continue;
        }

        std::vector<std::size_t> nonzIdx = tmp.b.nonzero_index();
        std::size_t dfsIdx = nonzIdx[dis(gen) % nonzIdx.size()];
        std::vector<ijkPair> ab = hgp_hrr_need(tmp, dfsIdx);
        for (const ijkPair &p : ab) {
            pairSet.insert(p);
        }
    }

    std::vector<ijkPair> pairVec(hrrPath.begin(), hrrPath.end());
    return pairVec;
}


void write_hgp_hrr_path(const std::vector<ijkPair> &path, std::size_t na, std::size_t nb) {
    std::string fileName = "hgp_hrr_path_" + std::to_string(na) + 
                                       "_" + std::to_string(nb) + 
                                       "_" + std::to_string(path.size()) + ".txt";
    std::ofstream ofs(fileName);
    assert(ofs.good());

    ofs << na << "    " << nb << std::endl;
    for (const ijkPair &p :  path) {
        ofs << p.to_string() << std::endl;
    }
    ofs.close();
}



int main(int argc, char **argv) {
    
    if (argc != 3) {
        std::cerr << "error argumnets." << std::endl;
        return -1;
    }
    
    std::size_t nTest = 500000;
    
    std::size_t na = std::stoi(argv[1]);
    std::size_t nb = std::stoi(argv[2]);

    if (na > 8 || nb > 8 || nb > na) {
        std::cerr << "error argumnets." << std::endl;
        return -1;
    }

    for (std::size_t n = 0; n < nTest; ++n) {
        std::printf("na = %2zu    nb = %2zu    try case: %5zu / %5zu\n",
                                                        na, nb, n, nTest);
        std::vector<ijkPair> pathThisTest = hgp_hrr_path_gen_random(na, nb);

        if (n == 0 || pathThisTest.size() < hgpHrrPathMin.size()) {
            hgpHrrPathMin = pathThisTest;
        }
    }

    std::sort(hgpHrrPathMin.begin(), hgpHrrPathMin.end(), [](const ijkPair &pa, const ijkPair &pb){
        if (pa.b.cnt_zero() != pb.b.cnt_zero()) return pa.b.cnt_zero() > pb.b.cnt_zero();
        return pa < pb;
    });

    write_hgp_hrr_path(hgpHrrPathMin, na, nb);

    return 0;
}



