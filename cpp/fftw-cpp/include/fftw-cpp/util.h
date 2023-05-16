#pragma once

#include <fftw3.h>
#include <cstddef>
#include <complex>
#include <numeric>
#include <concepts>

namespace fftw {

    // not enum classes for conversion and because namespacing isn't required
    enum Direction {
        FORWARD = FFTW_FORWARD,
        BACKWARD = FFTW_BACKWARD,
    };

    enum Flags {
        ESTIMATE = FFTW_ESTIMATE,
        MEASURE = FFTW_MEASURE,
        PATIENT = FFTW_PATIENT,
    };

    using std::size_t;

    namespace detail {
        // TODO specialize for float, long double, __float128
        template<std::floating_point Real>
        struct fftw_types;

        template<>
        struct fftw_types<double> {
            using complex = fftw_complex;
            using plan = fftw_plan;
        };

        template<std::floating_point Real>
        using fftw_complex_t = typename fftw_types<Real>::complex;

        template<std::floating_point Real>
        using fftw_plan_t = typename fftw_types<Real>::plan;


    } // fftw::detail

    template<bool IsReal, class Real, class Complex>
    using underlying_element_type = std::conditional_t<IsReal, Real, detail::fftw_complex_t<Real>>;

}