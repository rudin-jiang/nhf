#include "cg_base/hash.hpp"
#include "cg_base/ijk.hpp"
#include "cg_base/ijkpair.hpp"
#include "cg_base/utility.hpp"
#include <algorithm>
#include <unordered_set>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <string>
#include <fstream>
#include <utility>
#include <vector>
#include <cassert>


std::size_t na = 0;
std::size_t nb = 0;
std::vector<ijkPair> hgpHrrPath;

void read_hgp_path(std::string fileName) {
    std::ifstream ifs(fileName);
    assert(ifs.good());
    ifs >> na >> nb;

    std::size_t ai, aj, ak, bi, bj, bk;
    while (ifs >> ai >> aj >> ak >> bi >> bj >> bk) {
        ijk a(ai, aj, ak);
        ijk b(bi, bj, bk);
        hgpHrrPath.push_back(ijkPair(a, b));
    }

    ifs.close();
}

std::string hgp_hrr_code_filename(std::size_t na, std::size_t nb) {
    assert(na <= 8);
    assert(nb <= 8);
    return "hgp_hrr_" + encode_number(na) + "_" + encode_number(nb) + ".cpp";
}

std::string hgp_hrr_function_name(std::size_t na, std::size_t nb) {
    assert(na <= 8);
    assert(nb <= 8);
    return "hgp_hrr_" + encode_number(na) + "_" + encode_number(nb) + "_";
}

std::string hgp_hrr_function_signature(std::size_t na, std::size_t nb) {
    assert(na <= 8);
    assert(nb <= 8);
    return "void " + hgp_hrr_function_name(na, nb)
        + "(std::vector<double> &" + array_name(na, nb)
        + ", const std::vector<double> &hrrInp" 
        + ", double x, double y, double z)";
}

void print_head(std::ofstream &ofs) {
    // header file
    ofs << "#include \"hgp/hgp_hrr.hpp\"" << std::endl;
    ofs << "#include <vector>" << std::endl;
    ofs << std::endl;

    ofs << "namespace nhfInt {" << std::endl;
    ofs << "namespace hgp {" << std::endl;
    ofs << std::endl;

    ofs << hgp_hrr_function_signature(na, nb) << std::endl;
    ofs << "{" << std::endl;
    // ofs << std::endl;
}

void print_tail(std::ofstream &ofs) {
    // ofs << std::endl;
    ofs << "}" << std::endl;
    ofs << std::endl;
    ofs << "} // namespace (hgp)" << std::endl;
    ofs << "} // namespace (nhfInt)" << std::endl;
}


int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "error argument." << std::endl;
    }

    read_hgp_path(argv[1]);
    std::string codeFileName = hgp_hrr_code_filename(na, nb);
    std::ofstream ofs(codeFileName);
    assert(ofs.good());

    // handle the case of nb=0
    if (nb == 0) {
        print_head(ofs);

        ofs << "    for (std::size_t i = 0; i < " << shell_size(na) << "; ++i) {" << std::endl;
        ofs << "        " << array_name(na, nb) << "[i] = " << "hrrInp[i];" << std::endl;
        ofs << "    }" << std::endl;

        print_tail(ofs);
        ofs.close();
        return 0;
    }


    print_head(ofs);

    std::size_t cntMax = 0;
    for (std::size_t n = na; n < na + nb; ++n) {
        cntMax += (n+1) * (n+2) / 2;
    }
    std::size_t cntMaxLen = number_length(cntMax);

    std::size_t cnt = 0;
    for (std::size_t n = na; n <= na+nb; ++n) {
        ofs << "    ";
        ofs << "const double *" + array_name(n, 0)
            << " = hrrInp.data() + " + format_number(cnt, cntMaxLen)
            << ";" << std::endl;
        cnt += (n+1) * (n+2) / 2;
    }
    ofs << std::endl;

    // sort path
    auto cmp = [](const ijkPair &pa, const ijkPair &pb) {
        std::size_t paSumB = pa.b.sum();
        std::size_t pbSumB = pb.b.sum();
        if (paSumB != pbSumB) return paSumB < pbSumB;
        return pa < pb;
    };

    std::sort(hgpHrrPath.begin(), hgpHrrPath.end(), cmp);

    // get index max
    std::size_t indexMax = 0;
    for (const ijkPair &p : hgpHrrPath) {
        indexMax = std::max(indexMax, p.index());
    }
    std::size_t indexLen = number_length(indexMax + 1);

    // array need
    using SizePair = std::pair<std::size_t , std::size_t>;
    struct size_pair_hash {
        std::size_t operator()(const SizePair &p) const {
            return hash_val(p.first, p.second);
        }
    };

    std::unordered_set<SizePair, size_pair_hash> arrayNeed;
    for (const ijkPair &p : hgpHrrPath) {
        arrayNeed.insert({p.a.sum(), p.b.sum()});
    }

    std::vector<SizePair> sortedArrayNeed(
        arrayNeed.begin(), arrayNeed.end()
    );

    auto cmp2 = [](const SizePair &a, const SizePair &b){
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    };

    std::sort(sortedArrayNeed.begin(), sortedArrayNeed.end(), cmp2);

    std::size_t totSize = 0;

    for (const SizePair &p : sortedArrayNeed) {
        if (p.first == na && p.second == nb || p.second == 0) {
            continue;
        }

        std::size_t arraySize = pair_size(p.first, p.second);

        totSize += arraySize;

        ofs << "    "
            << "double "
            << array_name(p.first, p.second)
            << "[" + format_number(arraySize, indexLen) + "]"
            << ";"
            << std::endl;
    }
    ofs << std::endl;

    // std::cout << totSize << std::endl;

    // generate code
    std::vector<std::string> xyz = {"x", "y", "z"};
    std::unordered_set<ijkPair, ijkPair_hash> hasFormated;
    for (const ijkPair &p : hgpHrrPath) {
        if (p.b.cnt_zero() == 3) {
            hasFormated.insert(p);
            continue;
        }

        bool canBuild = false;
        for (std::size_t pos = 0; pos < 3; ++pos) {
            if (p.b[pos] == 0) {
                continue;
            }

            std::vector<ijkPair> needPair = p.hgp_hrr_need(pos);

            if (hasFormated.find(needPair[0]) != hasFormated.end() &&
                hasFormated.find(needPair[1]) != hasFormated.end()) {
                ofs << "    "
                    << p.array_element(indexLen)
                    << " = "
                    << needPair[0].array_element(indexLen)
                    << " + "
                    << xyz[pos]
                    << " * "
                    << needPair[1].array_element(indexLen)
                    << ";";
                
                ofs << "    //"
                    << p.to_comment()
                    << " = "
                    << needPair[0].to_comment()
                    << " + "
                    << "AB" + xyz[pos]
                    << " * "
                    << needPair[1].to_comment()
                    << std::endl;

                
                hasFormated.insert(p);
                canBuild = true;
                break;
            }
        }

        // std::cout << p.to_string() << std::endl;

        assert(canBuild);
    }

    print_tail(ofs);

    return 0;
}