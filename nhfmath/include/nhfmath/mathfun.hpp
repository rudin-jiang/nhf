#pragma once

#include <cstddef>
#include <cstdint>

namespace nhfMath {

// TODO: comments

// Maximum legal value of first parameter of factorial.
// In the current version: FactorialMaxN = 20
extern const std::size_t FactorialMaxN;

// Maximum legal value of first parameter of combination.
// In the current version: CombinationMaxN = 50
extern const std::size_t CombinationMaxN;

// Minimum and Maximum values of the legal 
// parameters of the function semifactorial.
// In the current version: SemifactorialMinN = -1
//                         SemifactorialMaxN = 30
extern const std::int64_t SemifactorialMinN;
extern const std::int64_t SemifactorialMaxN;


// n!
// https://en.wikipedia.org/wiki/Factorial
double factorial(std::size_t n);

// C(n,k)
// https://en.wikipedia.org/wiki/Combination
double combination(std::size_t n, std::size_t k);

// n!!
// https://en.wikipedia.org/wiki/Double_factorial
double semifactorial(std::int64_t n);

} // namespace (nhfMath)