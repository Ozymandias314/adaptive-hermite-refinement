#pragma once

#include <fftw3.h>
#include <complex>
#include "include_mdspan.h"
#include "util.h"

#include "basic_buffer_impl.h"
#include "basic_buffer.h"
#include "basic_plan.h"

namespace fftw {
    /// \defgroup Convenience types
    /// @{
    using plan = basic_plan<double>;
    using buffer = basic_buffer<double>;
    /// @}
} // fftw