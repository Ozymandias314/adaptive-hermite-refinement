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

            // Phi & Nabla
            for_each_kxky([&](Dim ky, Dim kx) {
                phi_K(kx, ky) = nonlinear::phi(moments_K(N_E, kx, ky), kx, ky);
                nablaPerpAPar_K(kx, ky) = kPerp(kx, ky) * moments_K(A_PAR, kx, ky);
                temp[0](kx, ky) = de * de * nablaPerpAPar_K(kx, ky);
                temp[1](kx, ky) = std::sqrt(2) * moments_K(G_MIN, kx, ky) + moments_K(N_E, kx, ky);
            });

            // store results of nonlinear operators, as well as results of predictor step
            fftw::mdbuffer<3u> GM_K_Star{M, X, Y}, GM_Nonlinear{M, X, Y};

            // Compute N
            auto bracketPhiNE_K = fullBracket(phi_K, sliceXY(moments_K, N_E));
            auto bracketAParNablaPerpAPar_K = fullBracket(sliceXY(moments_K, A_PAR), nablaPerpAPar_K);

            // Compute A
            auto bracketPhiAPar_K = fullBracket(phi_K, sliceXY(moments_K, A_PAR));
            auto bracketPhiDeNablaPerpAPar_K = fullBracket(phi_K, temp[0]);
            auto bracketNeG2APar_K = fullBracket(temp[1], sliceXY(moments_K, A_PAR));

            // Compute G2
            auto bracketPhiG2_K = fullBracket(phi_K, sliceXY(moments_K, G_MIN));
            auto bracketAParG3_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, 3));

            // Compute G_{M-1}
            auto bracketPhiGLast_K = fullBracket(phi_K, sliceXY(moments_K, LAST));
            auto bracketAParGLast_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, LAST));
            for_each_kxky([&](Dim kx, Dim ky) {
                bracketAParGLast_K(kx, ky) += rhoS / de * std::sqrt(LAST) * moments_K(LAST - 1, kx, ky);
            });
            auto bracketTotalGLast_K = fullBracket(sliceXY(moments_K, A_PAR), bracketAParGLast_K);

            for_each_kxky([&](Dim kx, Dim ky) {
                GM_Nonlinear(N_E, kx, ky) = nonlinear::N(bracketPhiNE_K(kx, ky), bracketAParNablaPerpAPar_K(kx, ky));
                GM_K_Star(N_E, kx, ky) = exp_nu(kx, ky, v2, dt) * moments_K(N_E, kx, ky) +
                                         dt / 2.0 * (1 + exp_nu(kx, ky, v2, dt)) * GM_Nonlinear(N_E, kx, ky);

                GM_Nonlinear(A_PAR, kx, ky) = nonlinear::A(bracketPhiAPar_K(kx, ky),
                                                           bracketPhiDeNablaPerpAPar_K(kx, ky),
                                                           bracketNeG2APar_K(kx, ky), kx, ky);
                GM_K_Star(A_PAR, kx, ky) = exp_eta(kx, ky, eta2, dt) * moments_K(A_PAR, kx, ky) +
                                           dt / 2.0 * (1 + exp_eta(kx, ky, eta2, dt)) * GM_Nonlinear(A_PAR, kx, ky);

                GM_Nonlinear(G_MIN, kx, ky) = nonlinear::G2(bracketPhiG2_K(kx, ky), bracketAParG3_K(kx, ky),
                                                            bracketAParNablaPerpAPar_K(kx, ky));
                GM_K_Star(G_MIN, kx, ky) = exp_nu(kx, ky, v2, dt) * moments_K(G_MIN, kx, ky) +
                                           dt / 2.0 * (1 + exp_nu(kx, ky, v2, dt)) * GM_Nonlinear(G_MIN, kx, ky);


                GM_Nonlinear(LAST, kx, ky) = nonlinear::GLast(bracketPhiGLast_K(kx, ky), bracketTotalGLast_K(kx, ky));
                GM_K_Star(LAST, kx, ky) =
                        exp_gm(LAST, hyper_nuei, dt) * exp_nu(kx, ky, v2, dt) * moments_K(LAST, kx, ky) +
                        dt / 2.0 * (1 + exp_gm(LAST, hyper_nuei, dt) * exp_nu(kx, ky, v2, dt)) *
                        GM_Nonlinear(LAST, kx, ky);
            });

            for (Dim m = 3; m < LAST; ++m) {
                auto bracketPhiGM_K = fullBracket(phi_K, sliceXY(moments_K, m));
                auto bracketAParGMMinus_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, m - 1));
                auto bracketAParGMPlus_K = fullBracket(sliceXY(moments_K, A_PAR), sliceXY(moments_K, m + 1));

                for_each_kxky([&](Dim kx, Dim ky) {
                    GM_Nonlinear(m, kx, ky) = nonlinear::GM(m, bracketPhiGM_K(kx, ky), bracketAParGMMinus_K(kx, ky),
                                                            bracketAParGMPlus_K(kx, ky));
                    GM_K_Star(m, kx, ky) = exp_gm(m, hyper_nuei, dt) * moments_K(m, kx, ky) +
                                           dt / 2.0 * (1 + exp_gm(m, hyper_nuei, dt)) * GM_Nonlinear(m, kx, ky);
                });
            }

            // TODO remove (currently preventing DCE)
            for_each_mxy([&](Dim m, Dim kx, Dim ky) {
                moments_K(m, kx, ky) = GM_K_Star(m, kx, ky);
            });

            // corrector step

            // Phi, Nabla, and other prep for A bracket
            for_each_kxky([&](Dim ky, Dim kx) {
                phi_K_Star(kx, ky) = nonlinear::phi(GM_K_Star(N_E, kx, ky), kx, ky);
                nablaPerpAPar_K_Star(kx, ky) = kPerp(kx, ky) * GM_K_Star(A_PAR, kx, ky);
                temp[0](kx, ky) = de * de * nablaPerpAPar_K_Star(kx, ky);
                temp[1](kx, ky) = std::sqrt(2) * GM_K_Star(G_MIN, kx, ky) + GM_K_Star(N_E, kx, ky);
            });

            // First, compute A_par
            auto bracketPhiAPar_K_Star = fullBracket(phi_K_Star, sliceXY(GM_K_Star, A_PAR));
            auto bracketPhiDeNablaPerpAPar_K_Star = fullBracket(phi_K_Star, temp[0]);
            auto bracketNeG2APar_K_Star = fullBracket(temp[1], sliceXY(GM_K_Star, A_PAR));

            fftw::mdbuffer<2u> nonlinearA_K_Star{X, Y}, semiImplicitOperator{X, Y};
            for_each_kxky([&](Dim kx, Dim ky) {

                nonlinearA_K_Star(kx, ky) = nonlinear::A(bracketPhiAPar_K_Star(kx, ky),
                                                         bracketPhiDeNablaPerpAPar_K_Star(kx, ky),
                                                         bracketNeG2APar_K_Star(kx, ky), kx, ky);

                semiImplicitOperator(kx, ky) = nonlinear::semiImplicitOp(dt, bPerpMax, aa0, kx, ky);
            });
// TODO
//            for (int p = 0; p <= MaxP; ++p) {
//
//            }
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
        DxDy<fftw::mdbuffer<2u>> derOp1{X, Y}, derOp2{X, Y};
        derivatives(op1, derOp1);
        derivatives(op2, derOp2);

        return halfBracket(derOp1, derOp2);
    }

    void Naive::derivatives(const Naive::ViewXY &op, Naive::DxDy<Naive::ViewXY> output) {
        DxDy<fftw::mdbuffer<2u>> Der_K{X, Y};
        prepareDXY_PH(op, Der_K.DX, Der_K.DY);
        fftInv(Der_K.DX.to_mdspan(), output.DY);
        fftInv(Der_K.DY.to_mdspan(), output.DY);
    }

    fftw::mdbuffer<2u> Naive::halfBracket(Naive::DxDy<Naive::ViewXY> derOp1, Naive::DxDy<Naive::ViewXY> derOp2) {
        fftw::mdbuffer<2u> br{X, Y}, br_K{X, Y};
        bracket(derOp1, derOp2, br);
        fft(br, br_K);

        return br_K;
    }
}