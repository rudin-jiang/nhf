#pragma once

#include <cstddef>
#include <vector>

namespace nhfInt {
namespace boysFun {


// Maximum legal value of first parameter of boysfun.
// In the current version: BoysFunMaxN = 64
extern const std::size_t BoysFunMaxN;

// The point for switching algorithms.
// when x <= BoysFunSwitch, evaluate using Taylor expansion
// when x >  BoysFunSwitch, evaluate using Asymptotic expression
// In the current version: BoysFunSwitch = 200.0
// (This value is exposed to facilitate testing and benchmark.)
extern const double BoysFunSwitch;

// Boys Function, see Helgaker et al. 
// Molecular electronic-structure theory. 
// John Wiley & Sons, 2013. Section 9.8
double boysfun(std::size_t n, double x);

// TODO:
// simultaneously generate a batch of boysfun
void boysfun(std::vector<double> &result, std::size_t n, double x);


} // namespace (boysFun)
} // namespace (nhfInt)