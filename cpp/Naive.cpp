#include "Naive.h"
#include "equillibrium.h"

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
    Naive::Naive(std::ostream &out, Dim M, Dim X, Dim Y) : HermiteRunner(out), M(M), X(X), Y(Y) {}

    void Naive::init(Dim N_) {
        N = N_;

        // Currently assuming X==Y for simplicity, but the code is written generally for the most part.
        assert(X == Y);

        // Plan FFTs both ways
        fft = fftw::plan<2u>::dft(phi_K, phi_K_New, fftw::FORWARD, fftw::MEASURE);
        fftInv = fftw::plan<2u>::dft(phi_K_New, phi_K, fftw::BACKWARD, fftw::MEASURE);

        // Initialize equilibrium values
        auto [aParEq, phi] = equilibriumOT01(X, Y);
        fft(phi, phi_K);
        fft(aParEq, aParEq_K);

        // DEBUG
        print(aParEq);
        print(aParEq_K);
        // Transform moments into phase space
        for (int m = 0; m < M; ++m) {
            if (m == A_PAR) continue;
            for_each_xy([&](Dim x, Dim y) {
                moments_K(m, x, y) = 0;
            });
        }

        // aPar and uekPar
        for_each_kxky([&](Dim ky, Dim kx) {
            aParEq_K(kx, ky) /= std::sqrt(KX*KY);
            moments_K(A_PAR, kx, ky) = aParEq_K(kx, ky);
            ueKPar_K(kx, ky) = -kPerp(kx, ky) * moments_K(A_PAR, kx, ky);
        });
        print(aParEq_K);
    }

    void Naive::run() {
        // store all derivatives
        DxDy<fftw::mdbuffer<3u>> dGM{M, X, Y};
        DxDy<fftw::mdbuffer<2u>> dPhi{X, Y}, dUEKPar{X, Y};

        bool divergent = false, repeat = false, noInc = false;
        int divergentCount = 0, repeatCount = 0;
        Real dt{0};
        HyperCoefficients hyper{};

        for (int t = 0; t < N; ++t) {
            // predictor step
            derivatives(phi_K, dPhi);
            derivatives(ueKPar_K, dPhi);
            for (int m = 0; m < M; ++m) {
                derivatives(sliceXY(moments_K, m), sliceXY(dGM, m));
            }

            // DEBUG
            print(sliceXY(moments_K, A_PAR));

            if (repeat or divergent) {
                t--; // repeat previous timestep
            } else {
                dt = getTimestep(dPhi, sliceXY(dGM, N_E), sliceXY(dGM, A_PAR));
                hyper = HyperCoefficients::calculate(dt, KX, KY, M);
            }

            // DEBUG
            std::cout << "dt: " << dt << std::endl;

            // store results of nonlinear operators, as well as results of predictor step
            fftw::mdbuffer<3u> GM_K_Star{M, X, Y}, GM_Nonlinear_K{M, X, Y};

            // Compute N
            auto bracketPhiNE_K = halfBracket(dPhi, sliceXY(dGM, N_E));
            auto bracketAParUEKPar_K = halfBracket(sliceXY(dGM, A_PAR), dUEKPar);

            // Compute A
            DxDy<fftw::mdbuffer<2u>> dPhiNeG2{X, Y};
            for_each_xy([&](Dim x, Dim y) {
                dPhiNeG2.DX(x, y) =
                        dPhi.DX(x, y) - rhoS * rhoS * (std::sqrt(2) * dGM.DX(G_MIN, x, y) + dGM.DX(N_E, x, y));
                dPhiNeG2.DY(x, y) =
                        dPhi.DY(x, y) - rhoS * rhoS * (std::sqrt(2) * dGM.DY(G_MIN, x, y) + dGM.DY(N_E, x, y));
            });

            auto bracketAParPhiG2Ne_K = halfBracket(sliceXY(dGM, A_PAR), dPhiNeG2);
            auto bracketPhiDeUEKPar_K = halfBracket(dPhi, dUEKPar);

            // Compute G2
            auto bracketPhiG2_K = halfBracket(dPhi, sliceXY(dGM, G_MIN));
            auto bracketAParG3_K = halfBracket(sliceXY(dGM, A_PAR), sliceXY(dGM, G_MIN + 1));

            // Compute G_{M-1}
            auto bracketPhiGLast_K = halfBracket(dPhi, sliceXY(dGM, LAST));
            auto bracketAParGLast_K = halfBracket(sliceXY(dGM, A_PAR), sliceXY(dGM, LAST));
            for_each_kxky([&](Dim kx, Dim ky) {
                bracketAParGLast_K(kx, ky) *= nonlinear::GLastBracketFactor(M, kx, ky, hyper);
                bracketAParGLast_K(kx, ky) += rhoS / de * std::sqrt(LAST) * moments_K(LAST - 1, kx, ky);
                // TODO Viriato adds this after the derivative
            });

            DxDy<fftw::mdbuffer<2u>> dBrLast{X, Y};
            derivatives(bracketAParGLast_K, dBrLast);
            auto bracketTotalGLast_K = halfBracket(sliceXY(dGM, A_PAR), dBrLast);

            for_each_kxky([&](Dim kx, Dim ky) {
                GM_Nonlinear_K(N_E, kx, ky) = nonlinear::N(bracketPhiNE_K(kx, ky), bracketAParUEKPar_K(kx, ky));
                GM_K_Star(N_E, kx, ky) = exp_nu(kx, ky, hyper.nu_2, dt) * moments_K(N_E, kx, ky) +
                                         dt / 2.0 * (1 + exp_nu(kx, ky, hyper.nu_2, dt)) * GM_Nonlinear_K(N_E, kx, ky);

                GM_Nonlinear_K(A_PAR, kx, ky) = nonlinear::A(bracketAParPhiG2Ne_K(kx, ky),
                                                             bracketPhiDeUEKPar_K(kx, ky), kx, ky);
                GM_K_Star(A_PAR, kx, ky) = exp_eta(kx, ky, hyper.eta2, dt) * moments_K(A_PAR, kx, ky) +
                                           dt / 2.0 * (1 + exp_eta(kx, ky, hyper.eta2, dt)) *
                                           GM_Nonlinear_K(A_PAR, kx, ky) +
                                           (1.0 - exp_eta(kx, ky, hyper.eta2, dt)) * aParEq_K(kx, ky);

                GM_Nonlinear_K(G_MIN, kx, ky) = nonlinear::G2(bracketPhiG2_K(kx, ky), bracketAParG3_K(kx, ky),
                                                              bracketAParUEKPar_K(kx, ky));
                GM_K_Star(G_MIN, kx, ky) = exp_nu(kx, ky, hyper.nu_2, dt) * moments_K(G_MIN, kx, ky) +
                                           dt / 2.0 * (1 + exp_nu(kx, ky, hyper.nu_2, dt)) *
                                           GM_Nonlinear_K(G_MIN, kx, ky);


                GM_Nonlinear_K(LAST, kx, ky) = nonlinear::GLast(bracketPhiGLast_K(kx, ky), bracketTotalGLast_K(kx, ky));
                GM_K_Star(LAST, kx, ky) =
                        exp_gm(LAST, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_g, dt) * moments_K(LAST, kx, ky) +
                        dt / 2.0 * (1 + exp_gm(LAST, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_g, dt)) *
                        GM_Nonlinear_K(LAST, kx, ky);
            });

            DxDy<fftw::mdbuffer<2u>> dGMinusPlus{X, Y};
            for (Dim m = 3; m < LAST; ++m) {
                for_each_xy([&](Dim x, Dim y) {
                    dGMinusPlus.DX(x, y) = std::sqrt(m) * dGM.DX(m - 1, x, y) + std::sqrt(m + 1) * dGM.DX(m + 1, x, y);
                    dGMinusPlus.DY(x, y) = std::sqrt(m) * dGM.DY(m - 1, x, y) + std::sqrt(m + 1) * dGM.DY(m + 1, x, y);
                });

                auto bracketAParGMMinusPlus_K = halfBracket(sliceXY(dGM, A_PAR), dGMinusPlus);
                auto bracketPhiGM_K = halfBracket(dPhi, sliceXY(dGM, m));

                for_each_kxky([&](Dim kx, Dim ky) {
                    GM_Nonlinear_K(m, kx, ky) = nonlinear::GM(m, bracketPhiGM_K(kx, ky),
                                                              bracketAParGMMinusPlus_K(kx, ky));
                    GM_K_Star(m, kx, ky) =
                            exp_gm(m, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_g, dt) * moments_K(m, kx, ky) +
                            dt / 2.0 * (1 + exp_gm(m, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_g, dt)) *
                            GM_Nonlinear_K(m, kx, ky);
                });
            }

            // corrector step

            // Phi, Nabla, and other prep for A bracket
            for_each_kxky([&](Dim ky, Dim kx) {
                phi_K_New(kx, ky) = nonlinear::phi(GM_K_Star(N_E, kx, ky), kx, ky);
                ueKPar_K_New(kx, ky) = -kPerp(kx, ky) * GM_K_Star(A_PAR, kx, ky);
            });


            DxDy<fftw::mdbuffer<2u>> dPhi_Loop{X, Y}, dUEKPar_Loop{X, Y};
            DxDy<fftw::mdbuffer<3u>> dGM_Loop{M, X, Y};
            derivatives(phi_K_New, dPhi_Loop);
            derivatives(ueKPar_K_New, dUEKPar_Loop);

            for (int m = 0; m < M; ++m) {
                // TODO(OPT) not necessary if we bail (only up to G_MIN)
                derivatives(sliceXY(GM_K_Star, m), sliceXY(dGM_Loop, m));
            }

            // Corrector loop
            // TODO confirm that only m derivatives are needed at a time
            //  (if not, can always store one in a temporary buffer)

            fftw::mdbuffer<2u> guessAPar_K{X, Y}, semiImplicitOperator{X, Y};
            for_each_kxky([&](Dim kx, Dim ky) {
                guessAPar_K(kx, ky) = moments_K(A_PAR, kx, ky);
                semiImplicitOperator(kx, ky) = nonlinear::semiImplicitOp(dt, bPerpMax, aa0, kx, ky);
            });

            Real old_error = 0, relative_error = 0;

            for (int p = 0; p <= MaxP; ++p) {
                auto DerivateNewMoment = [&](Dim m) {
                    derivatives(sliceXY(momentsNew_K, m), sliceXY(dGM_Loop, m));
                };

                // First, compute A_par
                DxDy<fftw::mdbuffer<2u>> dPhiNeG2_Loop{X, Y};
                for_each_xy([&](Dim x, Dim y) {
                    dPhiNeG2_Loop.DX(x, y) = dPhi_Loop.DX(x, y) -
                                             rhoS * rhoS *
                                             (std::sqrt(2) * dGM_Loop.DX(G_MIN, x, y) + dGM_Loop.DX(N_E, x, y));
                    dPhiNeG2_Loop.DY(x, y) = dPhi_Loop.DY(x, y) -
                                             rhoS * rhoS *
                                             (std::sqrt(2) * dGM_Loop.DY(G_MIN, x, y) + dGM_Loop.DY(N_E, x, y));
                });

                auto bracketAParPhiG2Ne_K_Loop = halfBracket(sliceXY(dGM_Loop, A_PAR), dPhiNeG2_Loop);
                auto bracketPhiDeUEKPar_K_Loop = halfBracket(dPhi_Loop, dUEKPar_Loop);

                /// f_pred from Viriato
                fftw::mdbuffer<3u> GM_Nonlinear_K_Loop{M, X, Y};
                Real sumAParRelError = 0;
                for_each_kxky([&](Dim kx, Dim ky) {
                    GM_Nonlinear_K_Loop(A_PAR, kx, ky) = nonlinear::A(bracketAParPhiG2Ne_K_Loop(kx, ky),
                                                                      bracketPhiDeUEKPar_K_Loop(kx, ky), kx, ky);
                    // TODO(OPT) reuse star
                    momentsNew_K(A_PAR, kx, ky) = 1.0 / (1.0 + semiImplicitOperator(kx, ky) / 4.0) *
                                                  (exp_eta(kx, ky, hyper.eta2, dt) * moments_K(A_PAR, kx, ky) +
                                                   dt / 2.0 * exp_eta(kx, ky, hyper.eta2, dt) *
                                                   GM_Nonlinear_K(A_PAR, kx, ky) +
                                                   dt / 2.0 * GM_Nonlinear_K_Loop(A_PAR, kx, ky) +
                                                   (1.0 - exp_eta(kx, ky, hyper.eta2, dt)) * aParEq_K(kx, ky) +
                                                   semiImplicitOperator(kx, ky) / 4.0 * guessAPar_K(kx, ky));
                    ueKPar_K_New(kx, ky) = -kPerp(kx, ky) * momentsNew_K(A_PAR, kx, ky);

                    sumAParRelError += std::norm(momentsNew_K(A_PAR, kx, ky) - moments_K(A_PAR, kx, ky));
                });

                old_error = relative_error;
                relative_error = 0;
                for_each_kxky([&](Dim kx, Dim ky) {
                    relative_error = std::max(relative_error, std::abs(
                            semiImplicitOperator(kx, ky) / 4.0 * (momentsNew_K(A_PAR, kx, ky) - guessAPar_K(kx, ky))) /
                                                              std::sqrt(sumAParRelError / (Real(KX) * Real(KY))));
                });

                std::cout << "rel err:" << relative_error << std::endl;
                // TODO bail if relative error is large

                DerivateNewMoment(A_PAR);
                auto bracketPhiNE_K_Loop = halfBracket(dPhi_Loop, sliceXY(dGM_Loop, N_E));
                auto bracketAParUEKPar_K_Loop = halfBracket(sliceXY(dGM_Loop, A_PAR), dUEKPar_Loop);

                for_each_kxky([&](Dim kx, Dim ky) {
                    GM_Nonlinear_K_Loop(N_E, kx, ky) = nonlinear::N(bracketPhiNE_K_Loop(kx, ky),
                                                                    bracketAParUEKPar_K_Loop(kx, ky));
                    // TODO(OPT) reuse star
                    momentsNew_K(N_E, kx, ky) = exp_nu(kx, ky, hyper.nu_2, dt) * moments_K(N_E, kx, ky) +
                                                dt / 2.0 * (1 + exp_nu(kx, ky, hyper.nu_2, dt)) * GM_Nonlinear_K(N_E, kx, ky) +
                                                dt / 2.0 * GM_Nonlinear_K_Loop(N_E, kx, ky);

                    phi_K_New(kx, ky) = nonlinear::phi(momentsNew_K(N_E, kx, ky), kx, ky);
                });

                derivatives(phi_K_New, dPhi_Loop);
                DerivateNewMoment(N_E);

                // Compute G2
                auto bracketPhiG2_K_Loop = halfBracket(dPhi_Loop, sliceXY(dGM_Loop, G_MIN));
                auto bracketAParG3_K_Loop = halfBracket(sliceXY(dGM_Loop, A_PAR), sliceXY(dGM_Loop, G_MIN + 1));

                for_each_kxky([&](Dim kx, Dim ky) {
                    GM_Nonlinear_K_Loop(G_MIN, kx, ky) = nonlinear::G2(bracketPhiG2_K_Loop(kx, ky),
                                                                       bracketAParG3_K_Loop(kx, ky),
                                                                       bracketAParUEKPar_K_Loop(kx, ky));
                    // TODO(OPT) reuse star
                    momentsNew_K(G_MIN, kx, ky) = exp_nu(kx, ky, hyper.nu_2, dt) * moments_K(G_MIN, kx, ky) +
                                                  dt / 2.0 * (1 + exp_nu(kx, ky, hyper.nu_2, dt)) *
                                                  GM_Nonlinear_K(G_MIN, kx, ky) +
                                                  dt / 2.0 * GM_Nonlinear_K_Loop(G_MIN, kx, ky);
                });
                DerivateNewMoment(G_MIN);

                DxDy<fftw::mdbuffer<2u>> dGMinusPlus_Loop{X, Y};
                for (int m = G_MIN + 1; m < LAST; ++m) {
                    for_each_xy([&](Dim x, Dim y) {
                        dGMinusPlus_Loop.DX(x, y) =
                                std::sqrt(m) * dGM_Loop.DX(m - 1, x, y) + std::sqrt(m + 1) * dGM_Loop.DX(m + 1, x, y);
                        dGMinusPlus_Loop.DY(x, y) =
                                std::sqrt(m) * dGM_Loop.DY(m - 1, x, y) + std::sqrt(m + 1) * dGM_Loop.DY(m + 1, x, y);
                    });

                    auto bracketAParGMMinusPlus_K_Loop = halfBracket(sliceXY(dGM_Loop, A_PAR), dGMinusPlus_Loop);
                    auto bracketPhiGM_K_Loop = halfBracket(dPhi_Loop, sliceXY(dGM_Loop, m));

                    for_each_kxky([&](Dim kx, Dim ky) {
                        GM_Nonlinear_K_Loop(m, kx, ky) = nonlinear::GM(m, bracketPhiGM_K_Loop(kx, ky),
                                                                       bracketAParGMMinusPlus_K_Loop(kx, ky));
                        // TODO(OPT) reuse star
                        momentsNew_K(m, kx, ky) =
                                exp_gm(m, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_g, dt) * moments_K(m, kx, ky) +
                                dt / 2.0 * exp_gm(m, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_g, dt) *
                                GM_Nonlinear_K(m, kx, ky) +
                                dt / 2.0 * GM_Nonlinear_K_Loop(m, kx, ky);
                    });

                    DerivateNewMoment(m);
                }

                // Compute G_{M-1}
                auto bracketPhiGLast_K_Loop = halfBracket(dPhi, sliceXY(dGM, LAST));
                auto bracketAParGLast_K_Loop = halfBracket(sliceXY(dGM, A_PAR), sliceXY(dGM, LAST));
                for_each_kxky([&](Dim kx, Dim ky) {
                    bracketAParGLast_K_Loop(kx, ky) *= nonlinear::GLastBracketFactor(M, kx, ky, hyper);
                    bracketAParGLast_K_Loop(kx, ky) += rhoS / de * std::sqrt(LAST) * momentsNew_K(LAST - 1, kx, ky);
                    // TODO Viriato adds this after the derivative
                });

                DxDy<fftw::mdbuffer<2u>> dBrLast_Loop{X, Y};
                derivatives(bracketAParGLast_K_Loop, dBrLast_Loop);
                auto bracketTotalGLast_K_Loop = halfBracket(sliceXY(dGM_Loop, A_PAR), dBrLast_Loop);

                for_each_kxky([&](Dim kx, Dim ky) {
                    GM_Nonlinear_K_Loop(LAST, kx, ky) = nonlinear::GLast(bracketPhiGLast_K_Loop(kx, ky),
                                                                         bracketTotalGLast_K_Loop(kx, ky));
                    // TODO(OPT) reuse star
                    momentsNew_K(LAST, kx, ky) =
                            exp_gm(LAST, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_2, dt) * moments_K(LAST, kx, ky) +
                            dt / 2.0 * exp_gm(LAST, hyper.nu_ei, dt) * exp_nu(kx, ky, hyper.nu_2, dt) *
                            GM_Nonlinear_K(LAST, kx, ky) +
                            dt / 2.0 * GM_Nonlinear_K_Loop(LAST, kx, ky);
                });
                DerivateNewMoment(LAST);

                if (relative_error <= epsilon) break;

                if (p != 0 and relative_error / old_error > 1.0) {
                    divergent = true;
                    divergentCount++;
                    dt = low * dt;
                    break;
                }
                if (p == MaxP) {
                    // did not converge well enough
                    repeat = true;
                    repeatCount++;
                    dt = low * dt;
                    break;
                }

                for_each_kxky([&](Dim kx, Dim ky) {
                    guessAPar_K(kx, ky) = momentsNew_K(A_PAR, kx, ky);
                });
            }
            if (divergent) continue;
            if (repeat) {
                noInc = true;
                continue;
            }

            Real tempDt = getTimestep(dPhi_Loop, sliceXY(dGM_Loop, N_E), sliceXY(dGM_Loop, A_PAR));
            dt = updateTimestep(dt, tempDt, noInc, relative_error);
            noInc = false;

            // New values are now old. Old values will be overwritten in the next timestep.
            std::swap(moments_K, momentsNew_K);
            std::swap(phi_K, phi_K_New);
            std::swap(ueKPar_K, ueKPar_K_New);
        }

        std::cout << "repeat: " << repeatCount << std::endl <<
                  "divergent: " << divergentCount << std::endl;
    }

    Real Naive::updateTimestep(Real dt, Real tempDt, bool noInc, Real relative_error) const {
        Real inc_factor = noInc ? 1 : 1.08;

        if (relative_error < 0.8 * epsilon) dt *= inc_factor;
        dt = std::min(tempDt, dt);
        return dt;
    }

    mdarray<Real, dextents<Dim, 2u>> Naive::getFinalAPar() {
        fftw::mdbuffer<2u> buf{X, Y};
        fftInv(sliceXY(moments_K, A_PAR), buf.to_mdspan());

        mdarray<Real, dextents<Dim, 2u>> result{X, Y};
        for_each_xy([&](Dim x, Dim y) {
            result(x, y) = buf(x, y).real();
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

        for_each_kxky([&](Dim kx, Dim ky) {
            // TODO filter
        });

        return br_K;
    }
}