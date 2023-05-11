#pragma once

#include "constants.h"

namespace ahr {
    [[nodiscard]] inline Real kPerp(Dim kx, Dim ky) {
        auto dkx = Real(kx), dky = Real(ky);
        return dkx * dkx + dky * dky;
    }

    [[nodiscard]] inline Real exp_nu(Dim kx, Dim ky, Real niu2, Real dt) {
        return std::exp(-(niu * kPerp(kx, ky) + niu2 * std::pow(kPerp(kx, ky), hyper_order) * dt));
    }

    [[nodiscard]] inline Real exp_gm(Dim m, Real hyper_nuei, Real dt) {
        return 0.0; // TODO
    }

    [[nodiscard]] inline Real exp_eta(Dim kx, Dim ky, Real res2, Real dt) {
        return 0.0; // TODO
    }

    [[nodiscard]] inline Real Gamma0(Real x) {
        Real ax = std::abs(x);

        if (ax < 3.75) {
            Real y = std::pow(x / 3.75, 2);
            return 1.0 + y * (3.5156229 +
                              y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
        } else {
            Real y = 3.75 / ax;
            return (std::exp(ax) / std::sqrt(ax)) *
                   (0.39894228 + y * (0.1328592e-1 + y *
                                                     (0.225319e-2 + y * (-0.157565e-2 + y *
                                                                                        (0.916281e-2 +
                                                                                         y * (-0.2057706e-1 + y *
                                                                                                              (0.2635537e-1 +
                                                                                                               y *
                                                                                                               (-0.1647633e-1 +
                                                                                                                y *
                                                                                                                0.392377e-2))))))));
        }
    }

    namespace nonlinear {

        [[nodiscard]] inline Complex phi(Complex ne_K, Dim kx, Dim ky) {
            if (rhoI<smallRhoI) return - ne_K *kPerp(kx, ky); else return rhoI *rhoI* 0.5 / (Gamma0(kPerp(kx, ky)*rhoI *rhoI* 0.5) - 1.0) *ne_K;
        }

        [[nodiscard]] inline Complex semiImplicitOp(Real dt, Real bPerpMax, Real aa0, Dim kx, Dim ky) {
            if (rhoI <= smallRhoI) {
                return aa0 * aa0 * (1 + kPerp(kx, ky) * (3.0 / 4.0 * rhoI * rhoI + rhoS * rhoS)) *
                       std::pow(kPerp(kx, ky) * bPerpMax * dt, 2) / (1 + kPerp(kx, ky) * de * de);
            } else {
                return aa0 * aa0 * (3.0 * rhoS * rhoS - rhoI * rhoI / (Gamma0(0.5 * kPerp(kx, ky) * rhoI * rhoI) - 1)) *
                       std::pow(kPerp(kx, ky) * bPerpMax * dt, 2) / (1 + kPerp(kx, ky) * de * de);
            }
        }

        [[nodiscard]] inline Complex N(Complex bracketPhiNE_K, Complex bracketAParUEKPar_K) {
            return -bracketPhiNE_K + bracketAParUEKPar_K;
        }

        [[nodiscard]] inline Complex A(Complex bracketPhiAPar_K, Complex bracketPhiDeUEKPar_K,
                                       Complex bracketNeG2APar_K, Dim kx, Dim ky) {
            return (-bracketPhiAPar_K + bracketPhiDeUEKPar_K +
                    rhoS * rhoS * bracketNeG2APar_K) / (1 + de * de * kPerp(kx, ky));
        }

        [[nodiscard]] inline Complex G2(Complex bracketPhiG2_K, Complex bracketAParG3_K,
                                        Complex bracketAParUEKPar_K) {
            return -bracketPhiG2_K +
                   std::sqrt(3.0) * rhoS / de * bracketAParG3_K +
                   std::sqrt(2.0) * bracketAParUEKPar_K;
        }

        [[nodiscard]] inline Complex GM(Dim m, Complex bracketPhiGM_K, Complex bracketAParGMMinus_K,
                                        Complex bracketAParGMPlus_K) {
            return -bracketPhiGM_K +
                   std::sqrt(m) * rhoS / de * bracketAParGMMinus_K +
                   std::sqrt(m + 1) * rhoS / de * bracketAParGMPlus_K;
            // TODO verify, maybe this is m and m+1?
        }

        [[nodiscard]] inline Complex GLast(Complex bracketPhiGLast_K, Complex bracketTotalGLast_K) {
            return -bracketPhiGLast_K + bracketTotalGLast_K;
        }
    }
}