#pragma once
namespace fftw {

    template<class Real, class Complex = std::complex<Real>>
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

} // namespace fftw