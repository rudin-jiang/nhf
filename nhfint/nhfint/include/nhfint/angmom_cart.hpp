#pragma once

#include <cstddef>
#include <vector>
#include <string>

namespace nhfInt {

// The maximum angular momentum supported 
// by the Cartesian Gaussian basis function.
// In the current version: AngMomCartMaxL = 8
// Attention:
// The maximum angular momentum of the basis function
// you can use is only 6. The case of angular momentum
// of 8 is considered here in order to deal with the 
// derivative of integrals.
extern const std::size_t AngMomCartMaxL;

// AngMomCart:
// a class used to represent cartesian angular momentum.
// Always guarantee that l <= AngMomCartMaxL && l == lx + ly + lz.
class AngMomCart {
public:
    std::size_t l;
    std::size_t lx, ly, lz;

    AngMomCart();
    AngMomCart(std::size_t l, std::size_t lx, 
               std::size_t ly, std::size_t lz);           
    AngMomCart(std::string line);
    
    std::size_t index() const;
    std::string to_string() const;

    // returns the components of angular momentum
    // [0] --> lx    [1] --> ly    [2] --> lz
    std::size_t operator[](std::size_t pos) const;
};

// To respect the tradition, the basis functions in the same shell are
// arranged in reverse lexicographic order according to (lx, ly, lz).
bool operator< (const AngMomCart &a, const AngMomCart &b);
bool operator> (const AngMomCart &a, const AngMomCart &b);
bool operator<=(const AngMomCart &a, const AngMomCart &b);
bool operator>=(const AngMomCart &a, const AngMomCart &b);
bool operator==(const AngMomCart &a, const AngMomCart &b);
bool operator!=(const AngMomCart &a, const AngMomCart &b);

// Return all AngMomCart that: lx + ly + lz == l,
// The returned vector is sorted in ascending order.
std::vector<AngMomCart> generate_angmomcart_shell(std::size_t l);

// return component of the idx'th angular momentum in a shell.
std::size_t angmomcart_component_x(std::size_t l, std::size_t idx);
std::size_t angmomcart_component_y(std::size_t l, std::size_t idx);
std::size_t angmomcart_component_z(std::size_t l, std::size_t idx);

}   // namespace (nhfInt)