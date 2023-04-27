#pragma once
namespace fftw {

    template<class Real, class Complex = std::complex<Real>>
    class basic_buffer : protected basic_buffer_impl<Real, Complex> {
    protected:
        using super = basic_buffer_impl<Real, Complex>;
    public:

        /// This constructor allocates the buffer but doesn't initialize it
        explicit basic_buffer(size_t length);

        /// This constructor allocates the buffer and initializes all elements to value
        basic_buffer(size_t length, Complex value);

        /// Deleted copy constructor to disable copying
        basic_buffer(const basic_buffer &) = delete;

        basic_buffer(basic_buffer &&other) noexcept; ///< Move constructor
        basic_buffer &operator=(basic_buffer &&other) noexcept; ///< Move assignment
        void swap(basic_buffer &other) noexcept;

        /// \defgroup Container methods (for range-for and other stdlib compatibility)
        /// @{
        using super::data;

        Complex *begin() { return super::data(); } ///<
        const Complex *begin() const { return super::data(); } ///<
        Complex *end() { return super::data() + length; }  ///<
        const Complex *end() const { return super::data() + length; }

        [[nodiscard]] size_t size() const { return length; } ///<
        Complex &operator[](size_t index) { return super::data()[index]; } ///<
        const Complex &operator[](size_t index) const { return super::data()[index]; } ///<
        /// @}

        using super::unwrap;

    private:
        size_t length{0};
        using super::storage;
    };

    template<class Real, class Complex>
    basic_buffer<Real, Complex>::basic_buffer(std::size_t length) :length(length) {
        super::storage = reinterpret_cast<detail::fftw_complex_t<Real> *>(fftw_malloc(length * sizeof(Complex)));
    }

    template<class Real, class Complex>
    void basic_buffer<Real, Complex>::swap(basic_buffer &other) noexcept {
        std::swap(length, other.length);
        std::swap(super::storage, other.storage);
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


}