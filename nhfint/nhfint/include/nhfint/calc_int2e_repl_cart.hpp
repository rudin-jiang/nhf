#pragma once

#include "nhfint/basis.hpp"
#include "nhfint/controller.hpp"
#include "nhfint/eritensor.hpp"


namespace nhfInt {

// segmented basis set
EriTensorD0 int2e_repl_d0_cart(const SegmentBasisSet &bsSet, const EriController &ctrl);
EriTensorD1 int2e_repl_d1_cart(const SegmentBasisSet &bsSet, const EriController &ctrl);
EriTensorD2 int2e_repl_d2_cart(const SegmentBasisSet &bsSet, const EriController &ctrl);

// general basis set
EriTensorD0 int2e_repl_d0_cart(const GeneralBasisSet &bsSet, const EriController &ctrl);
EriTensorD1 int2e_repl_d1_cart(const GeneralBasisSet &bsSet, const EriController &ctrl);
EriTensorD2 int2e_repl_d2_cart(const GeneralBasisSet &bsSet, const EriController &ctrl);

} // namespace (nhfInt)