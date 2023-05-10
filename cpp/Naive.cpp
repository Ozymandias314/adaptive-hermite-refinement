#include "Naive.h"

#include <utility>

namespace {
    template<class T, class ...Args, size_t ...I>
    auto forward_to_array_impl(std::index_sequence<I...> indices, Args &&...args) {
        return std::array<T, indices.size()>{
                [&args...](size_t index) {
                    return T{std::forward<Args>(args)...};
                }(I)...
        };
    }


    template<class T, size_t Size, class ...Args>
    std::array<T, Size> forward_to_array(Args &&...args) {
        return forward_to_array_impl<T>(std::make_index_sequence<Size>(), std::forward<Args>(args)...);
    }
}

namespace ahr {
    Naive::Naive(std::ostream &out, Dim M, Dim X, Dim Y) : HermiteRunner(out), M(M), X(X), Y(Y), temp{
            forward_to_array<fftw::mdbuffer<2u>, N_TEMP_BUFFERS>(X, Y)} {}

    void Naive::init(mdspan<Real, dextents<Dim, 3u>> initialMoments, Dim N_, Real initialDT_) {
        initialDT = initialDT_;
        N = N_;
        assert(M == initialMoments.extent(0));
        assert(X == initialMoments.extent(1));
        assert(Y == initialMoments.extent(2));

        // Currently assuming X==Y for simplicity, but the code is written generally for the most part.
        assert(X == Y);

        // Plan FFTs both ways
        fft = fftw::plan<2u>::dft(sliceXY(moments, 0), sliceXY(moments_K, 0), fftw::FORWARD, fftw::MEASURE);
        fftInv = fftw::plan<2u>::dft(sliceXY(moments_K, 0), sliceXY(moments, 0), fftw::BACKWARD, fftw::MEASURE);

        // Initialize moments
        // TODO only ne and A_PAR
        for_each_mxy([&](Dim m, Dim x, Dim y) {
            moments(m, x, y) = initialMoments(m, x, y);
        });

        // Transform moments into phase space
        for (int m = 0; m < M; ++m) {
            fft(sliceXY(moments, m), sliceXY(moments_K, m));
        }
    }

    void Naive::run() {
        Real dt = initialDT;
        for (int t = 0; t < N; ++t) {
            // predictor step

            // Phi
            // TODO maybe compute PHI? Or maybe that happened in the previous timestep

            // Nabla
            for_each_xy([&](Dim kx, Dim ky) {
                nablaPerpAPar_K(kx, ky) = kPerp(kx, ky) * moments_K(A_PAR, kx, ky);
            });

            // store results of nonlinear operators
            fftw::mdbuffer<3u> GM_K_Star{M, X, Y};

            // Compute N
            auto bracketPhiNE_K = fullBracket(phi_K, sliceXY(moments_K, N_E));
            auto bracketAParNablaPerpAPar_K = fullBracket(sliceXY(moments_K, A_PAR), nablaPerpAPar_K);
//            NonlinearN(bracketPhiNE_K, bracketAParNablaPerpAPar_K, sliceXY(GM_K_Star, N_E));
            for_each_kxky([&](Dim kx, Dim ky) {
                GM_K_Star(N_E, kx, ky) = exp_nu(kx, ky, v2, dt) * moments_K(N_E, kx, ky) +
                                         dt / 2.0 * (1 + exp_nu(kx, ky, v2, dt)) *
                                         nonlinear::N(bracketPhiNE_K(kx, ky), bracketAParNablaPerpAPar_K(kx, ky));
            });

            // Compute A
            auto bracketPhiAPar_K = fullBracket(phi_K, sliceXY(moments_K, A_PAR));
            for_each_xy([&](Dim kx, Dim ky) {
                temp[0](kx, ky) = de * de * nablaPerpAPar_K(kx, ky);
            });

            auto bracketPhiDeNablaPerpAPar_K = fullBracket(phi_K, temp[0]);

            for_each_xy([&](Dim kx, Dim ky) {
                temp[1](kx, ky) = std::sqrt(2) * moments_K(2, kx, ky) + moments_K(N_E, kx, ky);
            });
            auto bracketNeG2APar_K = fullBracket(temp[1], sliceXY(moments_K, A_PAR));
            NonlinearA(bracketPhiAPar_K, bracketPhiDeNablaPerpAPar_K, bracketNeG2APar_K, sliceXY(GM_K_Star, A_PAR));

            // Compute G2
            auto bracketPhiG2_K = fullBracket(phi_K, sliceXY(moments_K, 2));
            auto bracketAParG3_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, 3));
            NonlinearG2(bracketPhiG2_K, bracketAParG3_K, bracketAParNablaPerpAPar_K, sliceXY(GM_K_Star, 2));

            for (int m = 3; m < M; ++m) {
                auto bracketPhiGM_K = fullBracket(phi_K, sliceXY(moments_K, m));
                auto bracketAParGMMinus_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, m - 1));
                auto bracketAParGMPlus_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, m + 1));
                NonlinearGM(m, bracketPhiGM_K, bracketAParGMMinus_K, bracketAParGMPlus_K, sliceXY(GM_K_Star, m));
            }

            // corrector step
        }
    }

    mdarray<Real, dextents<Dim, 3u>> Naive::getFinalValues() {
        for (int m = 0; m < M; ++m) {
            fftInv(sliceXY(moments_K, m), sliceXY(moments, m));
        }

        mdarray<Real, dextents<Dim, 3u>> result{M, X, Y};
        for_each_mxy([&](Dim m, Dim x, Dim y) {
            result(m, x, y) = moments(m, x, y).real();
        });

        return result;
    }

    [[nodiscard]] fftw::mdbuffer<2u> Naive::fullBracket(Naive::ViewXY op1, Naive::ViewXY op2) {
        auto bufs = forward_to_array<fftw::mdbuffer<2u>, 8u>(X, Y);
        prepareDXY_PH(op1, bufs[0], bufs[1]);
        fftInv(bufs[0], bufs[2]);
        fftInv(bufs[1], bufs[3]);
        prepareDXY_PH(op2, bufs[4], bufs[5]);
        fftInv(bufs[4], bufs[6]);
        fftInv(bufs[5], bufs[7]);

        bracket(bufs[2], bufs[3], bufs[6], bufs[7], bufs[1]);
        fft(bufs[1], bufs[0]);
        return std::move(bufs[0]);
    }
}