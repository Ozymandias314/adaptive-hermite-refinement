#pragma once
namespace fftw {

    template<class Real, class Complex = std::complex<Real>, bool IsReal = false>
    class basic_buffer;

    template<class Real, class Complex = std::complex<Real>>
    using basic_rbuffer = basic_buffer<Real, Complex, true>;

    template<typename Real, typename Extents,
            typename Complex = std::complex<Real>,
            typename Layout = MDSPAN::layout_right,
            bool IsReal = false>
    using basic_mdbuffer = MDSPAN::mdarray<std::conditional_t<IsReal, Real, Complex>, Extents, Layout, fftw::basic_buffer<Real, Complex, IsReal>>;

    template<typename Real, typename Extents,
            typename Complex = std::complex<Real>,
            typename Layout = MDSPAN::layout_right>
    using basic_rmdbuffer = basic_mdbuffer<Real, Extents, Complex, Layout, true>;

    template<bool IsReal, class Real, class Complex = std::complex<double>>
    using underlying_element_type = std::conditional_t<IsReal, Real, detail::fftw_complex_t<Real>>;


    template<class Real, class Complex, bool IsReal>
    class basic_buffer {
    public:
        using element_type = std::conditional_t<IsReal, Real, Complex>;
        using pointer = element_type *;
        using reference = element_type &;
        using const_pointer = const element_type *;
        using const_reference = const element_type &;

        using underlying_element_type = std::conditional_t<IsReal, Real, detail::fftw_complex_t<Real>>;

        /// This constructor allocates the buffer but doesn't initialize it
        explicit basic_buffer(size_t length);

        /// This constructor allocates the buffer and initializes all elements to value
        basic_buffer(size_t length, element_type value);

        /// Deleted copy constructor to disable copying
        basic_buffer(const basic_buffer &) = delete;

        basic_buffer(basic_buffer &&other) noexcept; ///< Move constructor
        basic_buffer &operator=(basic_buffer &&other) noexcept; ///< Move assignment
        void swap(basic_buffer &other) noexcept;

        /// \defgroup Container methods (for range-for and other stdlib compatibility)
        /// @{

        auto data() -> element_type *; ///<
        auto data() const -> const element_type *; ///<
        element_type *begin() { return data(); } ///<
        const element_type *begin() const { return data(); } ///<
        element_type *end() { return data() + length; }  ///<
        const element_type *end() const { return data() + length; }

        [[nodiscard]] size_t size() const { return length; } ///<
        element_type &operator[](size_t index) { return data()[index]; } ///<
        const element_type &operator[](size_t index) const { return data()[index]; } ///<
        /// @}

        [[nodiscard]] underlying_element_type *unwrap() { return storage; }

        [[nodiscard]] const underlying_element_type *unwrap() const { return storage; }

        ~basic_buffer();

    private:
        size_t length{0};
        underlying_element_type *storage{nullptr};
    };

    template<class Real, class Complex, bool IsReal>
    basic_buffer<Real, Complex, IsReal>::basic_buffer(std::size_t length) :length(length) {
        storage = reinterpret_cast<underlying_element_type *>(fftw_malloc(length * sizeof(underlying_element_type)));
    }

    template<class Real, class Complex, bool IsReal>
    void basic_buffer<Real, Complex, IsReal>::swap(basic_buffer &other) noexcept {
        std::swap(length, other.length);
        std::swap(storage, other.storage);
    }

    template<class Real, class Complex, bool IsReal>
    basic_buffer<Real, Complex, IsReal>::basic_buffer(size_t length, element_type value) : basic_buffer(length) {
        for (Complex &elem: *this) {
            elem = value;
        }
    }

    template<class Real, class Complex, bool IsReal>
    basic_buffer<Real, Complex, IsReal>::basic_buffer(basic_buffer &&other) noexcept {
        other.swap(*this);
    }

    template<class Real, class Complex, bool IsReal>
    basic_buffer<Real, Complex, IsReal> &basic_buffer<Real, Complex, IsReal>::operator=(basic_buffer &&other) noexcept {
        basic_buffer(std::move(other)).swap(*this);
        return *this;
    }


}

template<class Real, class Complex, bool IsReal>
fftw::basic_buffer<Real, Complex, IsReal>::~basic_buffer() {
    if (storage) fftw_free(storage);
    storage = nullptr;
}

template<class Real, class Complex, bool IsReal>
auto fftw::basic_buffer<Real, Complex, IsReal>::data() -> element_type * {
    return reinterpret_cast<element_type *>(storage);
}

template<class Real, class Complex, bool IsReal>
auto
fftw::basic_buffer<Real, Complex, IsReal>::data() const -> const element_type * {
    return reinterpret_cast<const element_type *>(storage);
}