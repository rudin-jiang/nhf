#include "nhf/molecule.hpp"
#include <cstddef>

namespace nhf {

std::size_t Molecule::n_atom() const {
    return atomList.size();
}


} // namespace (nhf)