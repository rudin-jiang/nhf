#include "nhfint/angmom_cart.hpp"
#include <gtest/gtest.h>
#include <string>
#include <cstddef>
#include <algorithm>
#include <vector>
#include <set>
#include <cstdio>
#include <iostream>

static std::vector<nhfInt::AngMomCart>
generate_angmomcart_shell_brute_force(std::size_t l) {
    assert(l <= nhfInt::AngMomCartMaxL);

    std::set<nhfInt::AngMomCart> amcSet;
    for (std::size_t lx = 0; lx <= l; ++lx) {
    for (std::size_t ly = 0; ly <= l; ++ly) {
    for (std::size_t lz = 0; lz <= l; ++lz) {
        if (lx + ly + lz == l) {
            amcSet.insert(nhfInt::AngMomCart(l, lx, ly, lz));
        }
    }}}

    std::size_t sizeCalc = (l+1)*(l+2) / 2;
    assert(amcSet.size() == sizeCalc);

    return std::vector<nhfInt::AngMomCart>(
        amcSet.begin(), amcSet.end()
    );
}

TEST(TestAngMomCart, TestConstructor) {
    // default constructor
    nhfInt::AngMomCart amc0;
    EXPECT_EQ(amc0.l , 0);
    EXPECT_EQ(amc0.lx, 0);
    EXPECT_EQ(amc0.ly, 0);
    EXPECT_EQ(amc0.lz, 0);

    // value constructor
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        for (std::size_t lx = 0; lx <= l; ++lx) {
        for (std::size_t ly = 0; ly <= l; ++ly) {
        for (std::size_t lz = 0; lz <= l; ++lz) {
            if (lx + ly + lz == l) {
                nhfInt::AngMomCart amc1(l, lx, ly, lz);
                EXPECT_EQ(amc1.l , l);
                EXPECT_EQ(amc1.lx, lx);
                EXPECT_EQ(amc1.ly, ly);
                EXPECT_EQ(amc1.lz, lz);
            }
        }}}
    }

    // string constructor
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        for (std::size_t lx = 0; lx <= l; ++lx) {
        for (std::size_t ly = 0; ly <= l; ++ly) {
        for (std::size_t lz = 0; lz <= l; ++lz) {
            if (lx + ly + lz == l) {
                std::string line = " "
                    + std::to_string(l ) + " "
                    + std::to_string(lx) + " "
                    + std::to_string(ly) + " "
                    + std::to_string(lz) + " ";
                nhfInt::AngMomCart amc2(line);
                EXPECT_EQ(amc2.l , l);
                EXPECT_EQ(amc2.lx, lx);
                EXPECT_EQ(amc2.ly, ly);
                EXPECT_EQ(amc2.lz, lz);
            }
        }}}
    }
}

TEST(TestAngMomCart, TestIndex) {
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        std::vector<nhfInt::AngMomCart> amcVec
            = generate_angmomcart_shell_brute_force(l);
        
        for(std::size_t i = 0; i < amcVec.size(); ++i) {
            EXPECT_EQ(amcVec[i].index(), i);
        }
    }
}

TEST(TestAngMomCart, TestToString) {
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        for (std::size_t lx = 0; lx <= l; ++lx) {
        for (std::size_t ly = 0; ly <= l; ++ly) {
        for (std::size_t lz = 0; lz <= l; ++lz) {
            if (lx + ly + lz == l) {
                nhfInt::AngMomCart amc(l, lx, ly, lz);
                std::string llStr = (l  > 9 ? "" : " ") + std::to_string(l );
                std::string lxStr = (lx > 9 ? "" : " ") + std::to_string(lx);
                std::string lyStr = (ly > 9 ? "" : " ") + std::to_string(ly);
                std::string lzStr = (lz > 9 ? "" : " ") + std::to_string(lz);
                
                std::string line = 
                    llStr + "  " + lxStr + "  " + lyStr + "  " + lzStr;

                EXPECT_TRUE(amc.to_string() == line);
            }
        }}}
    }
}

TEST(TestAngMomCart, TestOperatorSquareBracket) {
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        for (std::size_t lx = 0; lx <= l; ++lx) {
        for (std::size_t ly = 0; ly <= l; ++ly) {
        for (std::size_t lz = 0; lz <= l; ++lz) {
            if (lx + ly + lz == l) {
                nhfInt::AngMomCart amc(l, lx, ly, lz);
                EXPECT_EQ(amc[0], lx);
                EXPECT_EQ(amc[1], ly);
                EXPECT_EQ(amc[2], lz);
            }
        }}}
    }
}

TEST(TestAngMomCart, TestOperatorCmp) {
    for (std::size_t l0 = 0; l0 <= nhfInt::AngMomCartMaxL; ++l0) {
    for (std::size_t l1 = 0; l1 <= nhfInt::AngMomCartMaxL; ++l1) {
        std::vector<nhfInt::AngMomCart> amcVec0
            = generate_angmomcart_shell_brute_force(l0);
        std::vector<nhfInt::AngMomCart> amcVec1
            = generate_angmomcart_shell_brute_force(l1);

        for (const nhfInt::AngMomCart &amc0 : amcVec0) {
        for (const nhfInt::AngMomCart &amc1 : amcVec1) {
            bool amc0Small = false;
            bool amc1Small = false;

            if (amc0.l != amc1.l) {
                amc0Small = (amc0.l < amc1.l);
                amc1Small = (amc1.l < amc0.l);
            }
            else { // amc0.l == amc1.l
                
                if (amc0.lx != amc1.lx) {
                    amc0Small = (amc0.lx > amc1.lx);
                    amc1Small = (amc1.lx > amc0.lx);
                }
                else { // amc0.lx == amc1.lx

                    if (amc0.ly != amc1.ly) {
                        amc0Small = (amc0.ly > amc1.ly);
                        amc1Small = (amc1.ly > amc0.ly);
                    }
                    else { // amc0.ly == amc1.ly
                        amc0Small = (amc0.lz > amc1.lz);
                        amc1Small = (amc1.lz > amc0.lz);
                    }
                }
            }

            if (amc0 <  amc1)
                EXPECT_TRUE( amc0Small && !amc1Small);
            
            if (amc0 >  amc1)
                EXPECT_TRUE(!amc0Small &&  amc1Small);

            if (amc0 <= amc1)
                EXPECT_TRUE(   ( amc0Small && !amc1Small)
                            || (!amc0Small && !amc1Small));

            if (amc0 >= amc1)
                EXPECT_TRUE(   (!amc0Small &&  amc1Small)
                            || (!amc0Small && !amc1Small));

            if (amc0 == amc1)
                EXPECT_TRUE(!amc0Small && !amc1Small);

            if (amc0 != amc1)
                EXPECT_TRUE(   ( amc0Small && !amc1Small)
                            || (!amc0Small &&  amc1Small));
        }}
    }}
}


TEST(TestAngMomCart, TestGenerateShell) {
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        std::vector<nhfInt::AngMomCart> amcVec0
            = nhfInt::generate_angmomcart_shell(l);
        
        std::vector<nhfInt::AngMomCart> amcVec1
            = generate_angmomcart_shell_brute_force(l);
        
        EXPECT_EQ(amcVec0.size(), amcVec1.size());
        for (std::size_t i = 0; i < amcVec0.size(); ++i) {
            EXPECT_EQ(amcVec0[i].l , amcVec1[i].l );
            EXPECT_EQ(amcVec0[i].lx, amcVec1[i].lx);
            EXPECT_EQ(amcVec0[i].ly, amcVec1[i].ly);
            EXPECT_EQ(amcVec0[i].lz, amcVec1[i].lz);
        }
    }
}


TEST(TestAngMomCart,TestComponent) {
    for (std::size_t l = 0; l <= nhfInt::AngMomCartMaxL; ++l) {
        std::vector<nhfInt::AngMomCart> amcVec
            = generate_angmomcart_shell_brute_force(l);
        
        for (std::size_t i = 0; i < amcVec.size(); ++i) {
            EXPECT_EQ(nhfInt::angmomcart_component_x(l, i), amcVec[i].lx);
            EXPECT_EQ(nhfInt::angmomcart_component_y(l, i), amcVec[i].ly);
            EXPECT_EQ(nhfInt::angmomcart_component_z(l, i), amcVec[i].lz);
        }
    }
}

