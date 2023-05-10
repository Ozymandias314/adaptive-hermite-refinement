#pragma once

#include "typedefs.h"

namespace ahr {
    constexpr Real rhoS = 0.1; // TODO
    constexpr Real de = 1.2; // TODO

    /// Diffusion
    constexpr Real niu = 100;
    constexpr Dim hyper_order = 3;
    constexpr Dim hyper_order_g = 3;
    constexpr Dim hyper_morder = 3;


    /// Hyper-coefficients
    // TODO make this runtime
    constexpr Real v2 = 0.1;
}