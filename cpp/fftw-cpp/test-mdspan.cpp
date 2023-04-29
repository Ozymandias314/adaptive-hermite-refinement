#include "fftw-cpp/fftw-cpp.h"

#include <iostream>
#include <numbers>

namespace stdex = std::experimental;

void print(const fftw::mdbuffer<3u> &arr, double scale = 1.0) {
    for (int j = 0; j < arr.extent(1); ++j) {
        for (int i = 0; i < arr.extent(0); ++i) {
            for (int k = 0; k < arr.extent(2); ++k) {
                std::cout << arr(i, j, k) * scale << " ";
            }
            std::cout << "| ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
};


int main() {
    std::size_t N = 3, M = 3, P = 3;
    fftw::mdbuffer<3u> in{N, M, P}, out{N, M, P}, out2{N, M, P};

    auto slice = [&](fftw::mdbuffer<3u> &buf, size_t index) {
        return stdex::submdspan(buf.to_mdspan(), index, stdex::full_extent, stdex::full_extent);
    };

    auto p = fftw::plan<2u>::dft(slice(in, 0), slice(out, 0), fftw::FORWARD, fftw::Flags::ESTIMATE);
    auto pInv = fftw::plan<2u>::dft(slice(out, 0), slice(out2, 0), fftw::BACKWARD, fftw::Flags::ESTIMATE);

    for (int i = 0; i < in.extent(0); ++i) {
        for (int j = 0; j < in.extent(1); ++j) {
            for (int k = 0; k < in.extent(2); ++k) {
                in(i, j, k) = {
                        std::cos(2.0 * std::numbers::pi * double(1 + i * k + j / 2) / (2.0 * in.size())),
                        std::cos(2.0 * std::numbers::pi * double(1 + i * k * j + in.size() / 2) / (2.0 * in.size()))
                };
            }
        }
    }

    print(in);

    for (int i = 0; i < in.extent(0); ++i) {
        p(slice(in, i), slice(out, i));
        pInv(slice(out, i), slice(out2, in.extent(0) - i - 1));
    }
    print(out);
    print(out2);
    print(out2, 1.0 / (slice(in, 0).size()));


    for (int i = 0; i < in.extent(0); ++i) {
        for (int j = 0; j < in.extent(1); ++j) {
            for (int k = 0; k < in.extent(2); ++k) {
                if (std::abs(in(i, j, k) - out2(in.extent(0) - i - 1, j, k) / double(slice(in, 0).size())) > 1.0e-4) {
                    std::cout << "Difference at index " << i << "," << j << "," << k << std::endl;
                    return 1;
                }
            }
        }
    }


}