#pragma once

#include "typedefs.h"

namespace ahr {
    constexpr Real rhoI = 4.6; // TODO
    constexpr Real rhoS = 0.1; // TODO
    constexpr Real de = 1.2; // TODO
    constexpr Real smallRhoI = 26.3; // TODO


    /// Diffusion
    constexpr Real niu = 100;
    constexpr Dim hyper_order = 3;
    constexpr Dim hyper_order_g = 3;
    constexpr Dim hyper_morder = 3;



    /// Hyper-coefficients
    // TODO make this runtime
    constexpr Real v2 = 0.1;
    constexpr Real eta2 = 0.1;
    constexpr Real hyper_nuei = 0.2;

    // TODO
    constexpr Real bPerpMax = 0.3;
    constexpr Real aa0 = 0.3;
}