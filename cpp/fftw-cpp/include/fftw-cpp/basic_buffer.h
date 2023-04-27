#pragma once
namespace fftw {

    template<class Real, class Complex = std::complex<Real>>
    class basic_buffer;

    template<typename Real, typename Extents,
            typename Complex = std::complex<Real>,
            typename Layout = MDSPAN::layout_right>
    using basic_mdbuffer = MDSPAN::mdarray<Complex, Extents, Layout, fftw::basic_buffer<Real, Complex>>;

    template<class Real, class Complex>
    class basic_buffer {
    public:
        using element_type = Complex;
        using pointer = Complex *;
        using reference = Complex &;
        using const_pointer = const Complex *;
        using const_reference = const Complex &;

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

        Complex *data(); ///<
        const Complex *data() const; ///<
        Complex *begin() { return data(); } ///<
        const Complex *begin() const { return data(); } ///<
        Complex *end() { return data() + length; }  ///<
        const Complex *end() const { return data() + length; }

        [[nodiscard]] size_t size() const { return length; } ///<
        Complex &operator[](size_t index) { return data()[index]; } ///<
        const Complex &operator[](size_t index) const { return data()[index]; } ///<
        /// @}

        [[nodiscard]] detail::fftw_complex_t<Real> *unwrap() { return storage; }

        [[nodiscard]] const detail::fftw_complex_t<Real> *unwrap() const { return storage; }

        ~basic_buffer();

    private:
        size_t length{0};
        detail::fftw_complex_t<Real> *storage{nullptr};
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


}

template<class Real, class Complex>
fftw::basic_buffer<Real, Complex>::~basic_buffer() {
    if (storage) fftw_free(storage);
    storage = nullptr;
}

template<class Real, class Complex>
Complex *fftw::basic_buffer<Real, Complex>::data() {
    return reinterpret_cast<Complex *>(storage);
}

template<class Real, class Complex>
const Complex *fftw::basic_buffer<Real, Complex>::data() const {
    return reinterpret_cast<Complex *>(storage);
}