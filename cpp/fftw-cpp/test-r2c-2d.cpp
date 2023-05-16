#include "fftw-cpp/fftw-cpp.h"

#include <iostream>
#include <numbers>


void print(const auto &arr) {
    for (int j = 0; j < arr.extent(0); ++j) {
        for (int k = 0; k < arr.extent(1); ++k) {
            std::cout << arr(j, k) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
};

int main() {
    namespace stdex = std::experimental;

    std::size_t N = 4, M = 6, MK = 4;
    fftw::rmdbuffer<2u> in{N, M};
    fftw::mdbuffer<2u> out{N, MK};
    fftw::basic_mdbuffer<double, stdex::dextents<std::size_t, 2u>, std::complex<double>, stdex::layout_left, true> out2{
            M, N};

    auto p = fftw::plan_r2c<2u>::dft(in.to_mdspan(), out.to_mdspan(), fftw::Flags::ESTIMATE);
    auto pInv = fftw::plan_c2r<2u>::dft(out.to_mdspan(), out2.to_mdspan(), fftw::Flags::ESTIMATE);

    for (int j = 0; j < in.extent(0); ++j) {
        for (int k = 0; k < in.extent(1); ++k) {
            in(j, k) = std::cos(2.0 * std::numbers::pi * double((j+1) * k) / (2.0 * in.size()));
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