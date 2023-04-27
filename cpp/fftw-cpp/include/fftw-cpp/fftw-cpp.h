#pragma once

#include <fftw3.h>
#include <complex>
#include "include_mdspan.h"
#include "util.h"

#include "basic_buffer.h"
#include "basic_plan.h"

namespace fftw {

    using MDSPAN::extents;
    using MDSPAN::dextents;
    using MDSPAN::dynamic_extent;
    using MDSPAN::layout_right;
    using MDSPAN::layout_left;
    using MDSPAN::layout_stride;

    /// \defgroup Convenience types
    /// @{
    template <size_t D = 1u>
    using plan = basic_plan<D, double>;
    using buffer = basic_buffer<double>;

    template<size_t D,
            typename Layout = MDSPAN::layout_right>
    using mdbuffer = basic_mdbuffer<double, dextents<size_t, D>, std::complex<double>, Layout>;
    /// @}

} // fftw