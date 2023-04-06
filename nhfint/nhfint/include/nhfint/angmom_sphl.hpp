// TODO: change data type


#pragma once

#include "nhfint/angmom_cart.hpp"
#include <cstddef>
#include <vector>
#include <string>
#include <cstdint>


namespace nhfInt {

extern const std::size_t AngMomSphlMaxL;

class AngMomSphl {
public:
    std::size_t     l;
    std::int64_t    lz;

    AngMomSphl();
    AngMomSphl(std::size_t l, std::int64_t lz);

    std::size_t index() const;
    std::string to_string() const;
};

bool operator< (const AngMomSphl &a, const AngMomSphl &b);
bool operator> (const AngMomSphl &a, const AngMomSphl &b);
bool operator<=(const AngMomSphl &a, const AngMomSphl &b);
bool operator>=(const AngMomSphl &a, const AngMomSphl &b);
bool operator==(const AngMomSphl &a, const AngMomSphl &b);
bool operator!=(const AngMomSphl &a, const AngMomSphl &b);

std::vector<AngMomSphl> generate_angmomsphl_shell(std::size_t l);

}   // namespace (nhfInt)


