#pragma once

#include "nhfmath/vec3d.hpp"
#include <cstddef>
#include <string>


namespace nhf {

class Atom {
public:
    std::size_t         type;
    nhfMath::Vec3d      centre;
    std::string         basis;

    Atom();
    Atom(std::size_t type, nhfMath::Vec3d centre, std::string basis);
};

} // namespace (nhf)

