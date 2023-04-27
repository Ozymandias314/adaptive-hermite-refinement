#pragma once

// =================
// Check mdspan support
// =================

#ifdef FFTW_CPP_MDSPAN_NAMESPACE
namespace fftw {
    /// Namespace that contains mdspan/mdarray and helpers (extents, default_accessor, etc.).
    /// Given by user.
    namespace MDSPAN = FFTW_CPP_MDSPAN_NAMESPACE;
};
#elif __cpp_lib_mdspan >= 202207L
#include <mdspan>
#include <mdarray>

namespace fftw {
    /// Namespace that contains mdspan/mdarray and helpers (extents, default_accessor, etc.)
    namespace MDSPAN = ::std;
};
#elif __has_include(<experimental/mdspan>)

#include <experimental/mdspan>
#include <experimental/mdarray>

namespace fftw {
    /// Namespace that contains mdspan/mdarray and helpers (extents, default_accessor, etc.)
    namespace MDSPAN = ::std::experimental;
};
#else
#error mdspan support required
#endif