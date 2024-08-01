#include "fftw-cpp/fftw-cpp.h"

#include "cilk.h"
#include "timing.h"

#include <argparse/argparse.hpp>
#include <cilk/cilkscale.h>
#include <iostream>
#include <numbers>

int main(int argc, const char *argv[]) {
    namespace stdex = std::experimental;
    argparse::ArgumentParser arguments("ahr");

    arguments.add_argument("X").help("Size of FFT domain").scan<'i', int>().default_value(128);
    arguments.add_argument("Y").help("Size of FFT domain").scan<'i', int>().default_value(128);
    arguments.add_argument("-c", "--check").help("Check correctness").flag();

    try {
        arguments.parse_args(argc, argv);
    } catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << arguments;
        return 1;
    }

    const auto X = arguments.get<int>("X");
    const auto Y = arguments.get<int>("Y");
    const auto check = arguments.get<bool>("--check");

    std::cout << "Running transforms of size " << X << "x" << X << std::endl;

    fftw::mdbuffer<2u> in{X, Y}, out{X, Y}, in2{X, Y}, mid2{X, Y}, out2{X, Y}, in3{X, Y},
        mid3{X, Y}, out3{X, Y};

    fftw::plan<2u> plan_2d = fftw::plan<2u>::dft(in, out, fftw::FORWARD, fftw::MEASURE);

    auto sliceX = [&](auto &buf, size_t x) {
        return stdex::submdspan(buf.to_mdspan(), x, stdex::full_extent);
    };
    auto sliceY = [&](auto &buf, size_t y) {
        return stdex::submdspan(buf.to_mdspan(), stdex::full_extent, y);
    };

    fftw::plan<1u> plan_1d_x =
        fftw::plan<1u>::dft(sliceY(in2, 0), sliceY(mid2, 0), fftw::FORWARD, fftw::MEASURE);
    fftw::plan<1u> plan_1d_y =
        fftw::plan<1u>::dft(sliceX(mid2, 0), sliceX(out, 0), fftw::FORWARD, fftw::MEASURE);

    for (int x1 = 0; x1 < in.extent(0); ++x1) {
        for (int x2 = 0; x2 < in.extent(1); ++x2) {
            in(x1, x2) = {std::cos(2.0 * std::numbers::pi * (x1 * x2) / (X * X)),
                          std::sin(2.0 * std::numbers::pi * (x1 * x2) / (X * X))};
            in3(x1, x2) = in2(x1, x2) = in(x1, x2);
        }
    }

    std::cout << "Planning complete!" << std::endl;
    std::cout << "2D transform:" << std::endl;
    auto [_, wsp] = benchmark(
        std::cout, [&]() {},
        [&]() {
            plan_2d();

            // normalize
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    out(x, y) /= X * Y;
                }
            }
        });

    // Do equivalent 1D FFTs
    std::cout << "1D transform:" << std::endl;
    auto [_2, wsp2] = benchmark(
        std::cout, [&]() {},
        [&]() {
            for (int x = 0; x < X; x++) {
                plan_1d_y(sliceX(in2, x), sliceX(mid2, x));
            }

            for (int y = 0; y < Y; y++) {
                plan_1d_x(sliceY(mid2, y), sliceY(out2, y));
            }

            // normalize
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    out2(x, y) /= X * Y;
                }
            }
        });

    // Do equivalent 1D FFTs with fusion
    std::cout << "1D transform (fusion):" << std::endl;
    auto [_3, wsp3] = benchmark(
        std::cout, [&]() {},
        [&]() {
            for (int x = 0; x < X; x++) {
                plan_1d_y(sliceX(in3, x), sliceX(mid3, x));
                for (int y = 0; y < Y; y++) {
                    mid3(x, y) /= X * Y;
                }
            }

            for (int y = 0; y < Y; y++) {
                plan_1d_x(sliceY(mid3, y), sliceY(out3, y));
            }
        });

    if (check) {
        for (int x = 0; x < X; x++) {
            for (int y = 0; y < Y; y++) {
                if (std::abs(out(x, y) - out2(x, y)) > 1e-10) {
                    std::cerr << "Mismatch at (" << x << ", " << y << "): " << out(x, y) << " vs "
                              << out2(x, y) << std::endl;
                    return 1;
                }
                if (std::abs(out(x, y) - out3(x, y)) > 1e-10) {
                    std::cerr << "Mismatch at (" << x << ", " << y << "): " << out(x, y) << " vs "
                              << out3(x, y) << std::endl;
                    return 1;
                }
            }
        }
    }

    for (unsigned i = 0; i < N_TRIALS; i++) {
        wsp_dump(wsp[i], "fft");
    }


}