#pragma once

#include "nhf/atom.hpp"
#include <cstddef>
#include <vector>
#include <cstdint>

namespace nhf {

class Molecule {
public:
    std::int64_t charge;
    std::int64_t multip;

    std::vector<Atom>   atomList;

    Molecule();
    Molecule(const std::string &info);

    std::size_t n_atom() const;
};

} // namespace (nhf)

