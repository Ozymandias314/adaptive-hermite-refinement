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

    template <size_t D = 1u>
    using plan_r2c = basic_plan_r2c<D, double>;

    template <size_t D = 1u>
    using plan_c2r = basic_plan_c2r<D, double>;

    using buffer = basic_buffer<double>;
    using rbuffer = basic_rbuffer<double>;

    template<size_t D,
            typename Layout = MDSPAN::layout_right, bool IsReal = false>
    using mdbuffer = basic_mdbuffer<double, dextents<size_t, D>, std::complex<double>, Layout, IsReal>;

    template<size_t D,
            typename Layout = MDSPAN::layout_right>
    using rmdbuffer = mdbuffer<D, Layout, true>;
    /// @}

} // fftw