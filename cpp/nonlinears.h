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

    [[nodiscard]] inline Real exp_ng(Real ng, Real hyper_nuei, Real dt) {
        return 0.0; // TODO
    }

    [[nodiscard]] inline Real exp_eta(Dim kx, Dim ky, Real res2, Real dt) {
        return 0.0; // TODO
    }

    namespace nonlinear {
        [[nodiscard]] inline Complex N(Complex bracketPhiNE_K, Complex bracketAParNablaPerpAPar_K) {
            return -bracketPhiNE_K + bracketAParNablaPerpAPar_K;
        }

        [[nodiscard]] inline Complex A(Complex bracketPhiAPar_K, Complex bracketPhiDeNablaPerpAPar_K,
                                       Complex bracketNeG2APar_K, Dim kx, Dim ky) {
            return (-bracketPhiAPar_K + bracketPhiDeNablaPerpAPar_K +
                    rhoS * rhoS * bracketNeG2APar_K) / (1 + de * de * kPerp(kx, ky));
        }

        [[nodiscard]] inline Complex G2(Complex bracketPhiG2_K, Complex bracketAParG3_K,
                                        Complex bracketAParNablaPerpAPar_K) {
            return -bracketPhiG2_K +
                   std::sqrt(3.0) * rhoS / de * bracketAParG3_K +
                   std::sqrt(2.0) * bracketAParNablaPerpAPar_K;
        }

        [[nodiscard]] inline Complex GM(Dim m, Complex bracketPhiGM_K, Complex bracketAParGMMinus_K,
                                        Complex bracketAParGMPlus_K) {
            return -bracketPhiGM_K +
                   std::sqrt(m - 1) * rhoS / de * bracketAParGMMinus_K +
                   std::sqrt(m) * rhoS / de * bracketAParGMPlus_K;
        }
    }
}