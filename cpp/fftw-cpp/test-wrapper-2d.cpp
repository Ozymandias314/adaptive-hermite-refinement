#include "fftw-cpp/fftw-cpp.h"

#include <iostream>
#include <numbers>

void print(const fftw::mdbuffer<2u> &arr) {
    for (int j = 0; j < arr.extent(0); ++j) {
        for (int k = 0; k < arr.extent(1); ++k) {
            std::cout << arr(j, k) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
};

int main() {
    std::size_t N = 4, M = 6;
    fftw::mdbuffer<2u> in{N, M}, out{N, M}, out2{N, M};

    static_assert(fftw::appropriate_buffer<2u, double, std::complex<double>, decltype(in)>);

    auto p = fftw::plan<2u>::dft(in, out, fftw::FORWARD, fftw::Flags::ESTIMATE);
    auto pInv = fftw::plan<2u>::dft(out, out2, fftw::BACKWARD, fftw::Flags::ESTIMATE);

    for (int j = 0; j < in.extent(0); ++j) {
        for (int k = 0; k < in.extent(1); ++k) {
            in(j, k) = {
                    std::cos(2.0 * std::numbers::pi * double(j * k) / (2.0 * in.size())),
                    std::cos(2.0 * std::numbers::pi * double(j * k + in.size() / 2) / (2.0 * in.size()))
            };
        }
    }

    print(in);

    p();
    pInv();

    print(out);
    print(out2);

    for (int j = 0; j < out2.extent(0); ++j) {
        for (int k = 0; k < out2.extent(1); ++k) {
            std::cout << out2(j, k) / double(out2.size()) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}