#include "fftw-cpp/fftw-cpp.h"

#include <numbers>
#include <iostream>

void print(const fftw::buffer &arr) {
    for (auto elem : arr) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
};

int main() {
    std::size_t N = 8;
    fftw::buffer in(N), out(N), out2(N);

    auto p = fftw::plan<>::dft(in, out, fftw::FORWARD, fftw::Flags::ESTIMATE);
    auto pInv = fftw::plan<>::dft(out, out2, fftw::BACKWARD, fftw::Flags::ESTIMATE);

    for (int j = 0; j < N; ++j) {
        in[j] = {
                std::cos(2.0 * std::numbers::pi * double(j) / double(N)),
                std::cos(2.0 * std::numbers::pi * double(j+N/2) / N )
        };
    }

    p();
    pInv();

    print(in);
    print(out);
    print(out2);

    for (auto elem : out2) {
        std::cout << elem / double(N) << " ";
    }
}