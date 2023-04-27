#pragma once

#include "basic_buffer.h"
#include "util.h"

namespace fftw {

    /// This boolean checks that the buffer is appropriate for this type of plan.
    /// By default, it is false
    template<size_t D, class Real, class Complex, typename T>
    constexpr inline bool appropriate_buffer = false;

    // We allow basic_buffer for 1D transforms
    template<class Real, class Complex = std::complex<Real>>
    constexpr inline bool appropriate_buffer<1u, Real, Complex, basic_buffer<Real, Complex>> = true;

    // We always allow a multi-d buffer for the same number of dimensions
    template<size_t D, class Real, class Complex, typename Layout, typename ExtentsIndexType, ExtentsIndexType... I>
    constexpr inline bool appropriate_buffer<D, Real, Complex, basic_mdbuffer<Real, MDSPAN::extents<ExtentsIndexType, I...>, Complex, Layout>>
            = sizeof...(I) == D;

    template<size_t D, class Real, class Complex = std::complex<Real>, typename T, typename T2>
    concept appropriate_buffers = appropriate_buffer<D, Real, Complex, T> &&
                                  appropriate_buffer<D, Real, Complex, T2>;

    template<size_t D, class Real, class Complex = std::complex<Real>>
    class basic_plan {
    protected:

    public:
        using real_t = Real;
        using complex_t = Complex;

        basic_plan() noexcept = default;

        ~basic_plan();

        /// Deleted copy constructor to disable copying
        basic_plan(const basic_plan &) = delete;

        basic_plan(basic_plan &&other) noexcept; ///< Move constructor
        basic_plan &operator=(basic_plan &&other) noexcept; ///< Move assignment
        void swap(basic_plan &other) noexcept;

        /// Executes the plan with the buffers provided initially.
        void operator()();

        template<typename BufferIn, typename BufferOut>
        requires appropriate_buffers<D, Real, Complex, BufferIn, BufferOut>
        void operator()(BufferIn &in, BufferOut &out);

        /// Returns the underlying FFTW plan.
        detail::fftw_plan_t<Real> unwrap() { return plan; }

        /// \defgroup{planning utilities}
        template<typename BufferIn, typename BufferOut>
        requires appropriate_buffers<D, Real, Complex, BufferIn, BufferOut>
        static auto dft(BufferIn &in, BufferOut &out, Direction direction, Flags flags) -> basic_plan;

    private:
        detail::fftw_plan_t<Real> plan{nullptr};
    };


    template<size_t D, class Real, class Complex>
    basic_plan<D, Real, Complex>::~basic_plan() {
        if (plan != nullptr) fftw_destroy_plan(plan);
        plan = nullptr;
    }

    template<size_t D, class Real, class Complex>
    void basic_plan<D, Real, Complex>::swap(basic_plan &other) noexcept {
        std::swap(plan, other.plan);
    }

    template<size_t D, class Real, class Complex>
    basic_plan<D, Real, Complex>::basic_plan(basic_plan &&other) noexcept {
        other.swap(*this);
    }

    template<size_t D, class Real, class Complex>
    basic_plan<D, Real, Complex> &basic_plan<D, Real, Complex>::operator=(basic_plan &&other) noexcept {
        basic_plan(other).swap(*this);
        return *this;
    }

    template<size_t D, class Real, class Complex>
    void basic_plan<D, Real, Complex>::operator()() {
        fftw_execute(plan);
    }

    template<size_t D, class Real, class Complex>
    template<typename BufferIn, typename BufferOut>
    requires appropriate_buffers<D, Real, Complex, BufferIn, BufferOut>
    void basic_plan<D, Real, Complex>::operator()(BufferIn &in, BufferOut &out) {
        fftw_execute_dft(plan, in.unwrap(), out.unwrap());
    }

    /// used for a static_assert inside an else block of if constexpr
    template<class...> inline constexpr bool always_false = false;

    namespace detail {

        template<class Real, class Complex>
        detail::fftw_complex_t<Real> *unwrap(basic_buffer<Real, Complex> &buf) {
            return buf.unwrap();
        }

        // TODO this is just a fix until we get a proper mdbuffer
        template<typename Real, typename Extents,
                typename Complex, typename Layout>
        detail::fftw_complex_t<Real> *unwrap(basic_mdbuffer<Real, Extents, Complex, Layout> &buf) {
            // for now, this is what FFT expects
            // TODO also put this in other places
            static_assert(std::is_same_v<Layout, MDSPAN::layout_right>);
            return reinterpret_cast<detail::fftw_complex_t<Real> *>(buf.data());
        }

        template<size_t D, class Real> requires (D == 1u) &&std::same_as<Real, double>

        auto plan_dft(auto &in, auto &out, Direction direction, Flags flags) {
            return fftw_plan_dft_1d(in.size(), unwrap(in), unwrap(out), direction, flags);
        }

        template<size_t D, class Real> requires (D == 2u) &&std::same_as<Real, double>

        auto plan_dft(auto &in, auto &out, Direction direction, Flags flags) {
            // TODO for layout left this is different
            return fftw_plan_dft_2d(in.extent(0), in.extent(1), unwrap(in), unwrap(out), direction, flags);
        }
    }

    template<size_t D, class Real, class Complex>
    template<typename BufferIn, typename BufferOut>
    requires appropriate_buffers<D, Real, Complex, BufferIn, BufferOut>
    auto basic_plan<D, Real, Complex>::dft(BufferIn &in, BufferOut &out, Direction direction,
                                           Flags flags) -> basic_plan {
        if (in.size() != out.size()) throw std::invalid_argument("mismatched buffer sizes");
        if (direction != FORWARD and direction != BACKWARD)
            throw std::invalid_argument("invalid direction");

        basic_plan plan1;
        plan1.plan = detail::template plan_dft<D, Real>(in, out, direction, flags);
        return plan1;
    }

} // namespace fftw