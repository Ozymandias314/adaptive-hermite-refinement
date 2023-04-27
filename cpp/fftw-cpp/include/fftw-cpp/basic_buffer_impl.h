#pragma once

namespace fftw {
    template<class Real, class Complex>
    class basic_buffer_impl {
    protected:
        detail::fftw_complex_t<Real> *storage{nullptr};
    public:

        ~basic_buffer_impl();

        [[nodiscard]] detail::fftw_complex_t<Real> *unwrap() { return storage; }

        [[nodiscard]] const detail::fftw_complex_t<Real> *unwrap() const { return storage; }

        Complex *data();
        const Complex *data() const;
    };

    template<class Real, class Complex>
    basic_buffer_impl<Real, Complex>::~basic_buffer_impl() {
        if (storage) fftw_free(storage);
        storage = nullptr;
    }

    template<class Real, class Complex>
    Complex *basic_buffer_impl<Real, Complex>::data() {
        return reinterpret_cast<Complex *>(storage);
    }

    template<class Real, class Complex>
    const Complex *basic_buffer_impl<Real, Complex>::data() const {
        return reinterpret_cast<Complex *>(storage);
    }
} // namespace fftw