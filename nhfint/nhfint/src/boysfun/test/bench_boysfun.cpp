#include "boysfun/boysfun.hpp"
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>

class BM_BOYSFUN : public benchmark::Fixture {
public:
    virtual void SetUp(const ::benchmark::State&) override {
    }
};

BENCHMARK_F(BM_BOYSFUN, SmallX)(benchmark::State& state)
{
    std::string filename = "test_data/boysfun_test_data_small_x.dat";
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;
    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
        assert(xValueList[i] <= nhfInt::boysFun::BoysFunSwitch);
    }

    for (auto _ : state) {
        double notOpt = 0.0;
        for (const double x : xValueList) {
            for (std::size_t n = 0; n <= nhfInt::boysFun::BoysFunMaxN; ++n) {
                notOpt += nhfInt::boysFun::boysfun(n, x);
            }
        }
        benchmark::DoNotOptimize(notOpt);
    }
}


BENCHMARK_F(BM_BOYSFUN, LargeX)(benchmark::State& state)
{
    std::string filename = "test_data/boysfun_test_data_large_x.dat";
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;
    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
        assert(xValueList[i] > nhfInt::boysFun::BoysFunSwitch);
    }

    for (auto _ : state) {
        double notOpt = 0.0;
        for (const double x : xValueList) {
            for (std::size_t n = 0; n <= nhfInt::boysFun::BoysFunMaxN; ++n) {
                notOpt += nhfInt::boysFun::boysfun(n, x);
            }
        }
        benchmark::DoNotOptimize(notOpt);
    }
}

BENCHMARK_F(BM_BOYSFUN, SmallX_SmallN)(benchmark::State& state)
{
    std::string filename = "test_data/boysfun_test_data_small_x.dat";
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;
    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
        assert(xValueList[i] <= nhfInt::boysFun::BoysFunSwitch);
    }

    for (auto _ : state) {
        double notOpt = 0.0;
        for (const double x : xValueList) {
            for (std::size_t n = 0; n <= nhfInt::boysFun::BoysFunMaxN / 2; ++n) {
                notOpt += nhfInt::boysFun::boysfun(n, x);
            }
        }
        benchmark::DoNotOptimize(notOpt);
    }
}

BENCHMARK_F(BM_BOYSFUN, SmallX_LargeN)(benchmark::State& state)
{
    std::string filename = "test_data/boysfun_test_data_small_x.dat";
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;
    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
        assert(xValueList[i] <= nhfInt::boysFun::BoysFunSwitch);
    }

    for (auto _ : state) {
        double notOpt = 0.0;
        for (const double x : xValueList) {
            for (std::size_t n = nhfInt::boysFun::BoysFunMaxN / 2 + 1;
                                n <= nhfInt::boysFun::BoysFunMaxN; ++n) {
                notOpt += nhfInt::boysFun::boysfun(n, x);
            }
        }
        benchmark::DoNotOptimize(notOpt);
    }
}

BENCHMARK_F(BM_BOYSFUN, LargeX_SmallN)(benchmark::State& state)
{
    std::string filename = "test_data/boysfun_test_data_large_x.dat";
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;
    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
        assert(xValueList[i] > nhfInt::boysFun::BoysFunSwitch);
    }

    for (auto _ : state) {
        double notOpt = 0.0;
        for (const double x : xValueList) {
            for (std::size_t n = 0; n <= nhfInt::boysFun::BoysFunMaxN / 2; ++n) {
                notOpt += nhfInt::boysFun::boysfun(n, x);
            }
        }
        benchmark::DoNotOptimize(notOpt);
    }
}

BENCHMARK_F(BM_BOYSFUN, LargeX_LargeN)(benchmark::State& state)
{
    std::string filename = "test_data/boysfun_test_data_large_x.dat";
    std::ifstream ifs(filename);
    assert(ifs.good());

    std::size_t nTestCase = 0;
    ifs >> nTestCase;
    std::vector<double> xValueList(nTestCase);
    for (std::size_t i = 0; i < nTestCase; ++i) {
        ifs >> xValueList[i];
        assert(xValueList[i] > nhfInt::boysFun::BoysFunSwitch);
    }

    for (auto _ : state) {
        double notOpt = 0.0;
        for (const double x : xValueList) {
            for (std::size_t n = nhfInt::boysFun::BoysFunMaxN / 2 + 1;
                                n <= nhfInt::boysFun::BoysFunMaxN; ++n) {
                notOpt += nhfInt::boysFun::boysfun(n, x);
            }
        }
        benchmark::DoNotOptimize(notOpt);
    }
}