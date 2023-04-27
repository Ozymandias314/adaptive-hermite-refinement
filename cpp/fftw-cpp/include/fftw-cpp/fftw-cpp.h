#pragma once

#include <fftw3.h>
#include <complex>



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

namespace fftw {
    /// Namespace that contains mdspan/mdarray and helpers (extents, default_accessor, etc.)
    namespace MDSPAN = ::std;
};
#elif __has_include(<experimental/mdspan>)

#include <experimental/mdspan>

namespace fftw {
    /// Namespace that contains mdspan/mdarray and helpers (extents, default_accessor, etc.)
    namespace MDSPAN = ::std::experimental;
};
#endif

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

    template<class Real, class Complex = std::complex<Real>>
    class basic_plan;

    template<class Real, class Complex = std::complex<Real>>
    class basic_buffer;

    /// \defgroup Convenience types
    /// @{
    using plan = basic_plan<double>;
    using buffer = basic_buffer<double>;
    using planf = basic_plan<float>;
    using bufferf = basic_buffer<float>;
    using planl = basic_plan<long double>;
    using bufferl = basic_buffer<long double>;
    /// @}

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

    using std::size_t;

    template<class Real, class Complex>
    class basic_buffer {
    public:

        /// This constructor allocates the buffer but doesn't initialize it
        explicit basic_buffer(size_t length);

        /// This constructor allocates the buffer and initializes all elements to value
        basic_buffer(size_t length, Complex value);

        ~basic_buffer();

        /// Deleted copy constructor to disable copying
        basic_buffer(const basic_buffer &) = delete;

        basic_buffer(basic_buffer &&other) noexcept; ///< Move constructor
        basic_buffer &operator=(basic_buffer &&other) noexcept; ///< Move assignment
        void swap(basic_buffer &other) noexcept;

        /// \defgroup Container methods (for range-for and other stdlib compatibility)
        /// @{
        Complex *data(); ///<
        const Complex *data() const; ///<
        Complex *begin() { return data(); } ///<
        Complex *end() { return data() + length; }  ///<
        const Complex *begin() const { return data(); } ///<
        const Complex *end() const { return data() + length; }

        [[nodiscard]] size_t size() const { return length; } ///<
        Complex &operator[](size_t index) { return data()[index]; } ///<
        const Complex &operator[](size_t index) const { return data()[index]; } ///<
        /// @}

        detail::fftw_complex_t<Real> *unwrap() { return storage; } ///<
        const detail::fftw_complex_t<Real> *unwrap() const { return storage; }

    private:
        detail::fftw_complex_t<Real> *storage{nullptr};
        size_t length{0};
    };

    template<class Real, class Complex>
    class basic_plan {
    public:
        using real_t = Real;
        using complex_t = Complex;
        using buffer_t = basic_buffer<Real, Complex>;

        basic_plan() noexcept = default;

        ~basic_plan();

        /// Deleted copy constructor to disable copying
        basic_plan(const basic_plan &) = delete;

        basic_plan(basic_plan &&other) noexcept; ///< Move constructor
        basic_plan &operator=(basic_plan &&other) noexcept; ///< Move assignment
        void swap(basic_plan &other) noexcept;

        /// Executes the plan with the buffers provided initially.
        void operator()();

        void operator()(buffer_t &in, buffer_t &out);

        /// Returns the underlying FFTW plan.
        detail::fftw_plan_t<Real> unwrap() { return plan; }

        /// \defgroup{planning utilities}
        template<std::size_t Dimension = 1u>
        static basic_plan dft(buffer_t &in, buffer_t &out, Direction direction, Flags flags);

    private:
        detail::fftw_plan_t<Real> plan{nullptr};
    };

    template<class Real, class Complex>
    basic_buffer<Real, Complex>::basic_buffer(std::size_t length) :length(length) {
        storage = reinterpret_cast<detail::fftw_complex_t<Real> *>(fftw_malloc(length * sizeof(Complex)));
    }

    template<class Real, class Complex>
    void basic_buffer<Real, Complex>::swap(basic_buffer &other) noexcept {
        std::swap(length, other.length);
        std::swap(storage, other.storage);
    }

    template<class Real, class Complex>
    basic_buffer<Real, Complex>::basic_buffer(size_t length, Complex value) : basic_buffer(length) {
        for (Complex &elem: *this) {
            elem = value;
        }
    }

    template<class Real, class Complex>
    basic_buffer<Real, Complex>::basic_buffer(basic_buffer &&other) noexcept {
        other.swap(this);
    }

    template<class Real, class Complex>
    basic_buffer<Real, Complex> &basic_buffer<Real, Complex>::operator=(basic_buffer &&other) noexcept {
        basic_buffer(other).swap(*this);
        return *this;
    }

    template<class Real, class Complex>
    basic_buffer<Real, Complex>::~basic_buffer() {
        if (storage) fftw_free(storage);
        storage = nullptr;
    }

    template<class Real, class Complex>
    Complex *basic_buffer<Real, Complex>::data() {
        return reinterpret_cast<Complex *>(storage);
    }

    template<class Real, class Complex>
    const Complex *basic_buffer<Real, Complex>::data() const {
        return reinterpret_cast<Complex *>(storage);
    }

    template<class Real, class Complex>
    basic_plan<Real, Complex>::~basic_plan() {
        if (plan != nullptr) fftw_destroy_plan(plan);
        plan = nullptr;
    }

    template<class Real, class Complex>
    void basic_plan<Real, Complex>::swap(basic_plan &other) noexcept {
        std::swap(plan, other.plan);
    }

    template<class Real, class Complex>
    basic_plan<Real, Complex>::basic_plan(basic_plan &&other) noexcept {
        other.swap(*this);
    }

    template<class Real, class Complex>
    basic_plan<Real, Complex> &basic_plan<Real, Complex>::operator=(basic_plan &&other) noexcept {
        basic_plan(other).swap(*this);
        return *this;
    }

    template<class Real, class Complex>
    void basic_plan<Real, Complex>::operator()() {
        fftw_execute(plan);
    }

    template<class Real, class Complex>
    void basic_plan<Real, Complex>::operator()(buffer_t &in, buffer_t &out) {
        fftw_execute_dft(plan, in.unwrap(), out.unwrap());
    }

    /// used for a static_assert inside an else block of if constexpr
    template<class...> inline constexpr bool always_false = false;

    template<class Real, class Complex>
    template<std::size_t Dimension>
    basic_plan<Real, Complex>
    basic_plan<Real, Complex>::dft(buffer_t &in, buffer_t &out, Direction direction, Flags flags) {
        basic_plan plan1;
        if (in.size() != out.size()) throw std::invalid_argument("mismatched buffer sizes");
        if (direction != FORWARD and direction != BACKWARD) throw std::invalid_argument("invalid direction");
        if constexpr (Dimension == 1u) {
            plan1.plan = fftw_plan_dft_1d(in.size(), in.unwrap(), out.unwrap(), direction, flags);
        } else {
            // needs to depend on some template param to not get instantiated.
            static_assert(always_false<Real> && "unimplemented dimension");
        }

        return plan1;
    }
} // fftw