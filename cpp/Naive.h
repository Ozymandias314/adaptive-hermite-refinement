#pragma once

#include "HermiteRunner.h"
#include "constants.h"
#include "nonlinears.h"

#include <fftw-cpp/fftw-cpp.h>
#include <type_traits>

namespace ahr {
    namespace stdex = std::experimental;

    class Naive : public ahr::HermiteRunner {
    public:
        /**
         * Construct a new simulation using the Naive approach.
         * @param out The stream to which we're logging.
         * @param M The number of moments. Note that moments are 0-indexed, meaning that the highest moment is M-1.
         * @param X The size of the X domain.
         * @param Y The size of the Y domain.
         */
        Naive(std::ostream &out, Dim M, Dim X, Dim Y);

        void init(Dim N_) override;

        void run() override;

        mdarray<Real, dextents<Dim, 2u>> getFinalAPar() override;

        using ViewXY = stdex::mdspan<Complex, stdex::dextents<Dim, 2u>>;
    private:
        Dim const M, X, Y, KX{X}, KY{Y};
        Dim N{};
        fftw::plan<2u> fft{}, fftInv{};
        Real bPerpMax{0};

        static constexpr Dim N_E = 0;
        static constexpr Dim A_PAR = 1;
        static constexpr Dim G_MIN = 2;
        const Dim LAST = M - 1;

        template<class Buffer>
        struct DxDy {
            using buffer_t = Buffer;
            Buffer DX, DY;

            template<typename ...Args>
            requires std::constructible_from<Buffer, Args...>
            explicit DxDy(Args &&...args) :
                    DX{std::forward<Args>(args)...},
                    DY{std::forward<Args>(args)...} {}

            DxDy(Buffer dx, Buffer dy) : DX(dx), DY(dy) {}

            template<class U>
            requires std::convertible_to<Buffer, U>
            operator DxDy<U>() { // NOLINT(google-explicit-constructor)
                return {U(DX), U(DY)};
            }
        };



        /// \defgroup Buffers for all the physical quantities used.
        /// Names ending in PH mean the values are in phase space.
        /// @{

        /// g_m: moment values for moments $m \in [0,M-1]$.
        /// The following values are special moments:
        /// - m=0: n_e (charge density)
        /// - m=1: A∥ (or Apar, parallel velocity)
        /// Additionally, moments is used only at the start and end.
        /// During the simulation, we use moments
        /// TODO maybe instead of these enormous amounts of memory, we could reuse (parallelism might suffer)
        fftw::mdbuffer<3u> moments_K{M, X, Y}, momentsNew_K{M, X, Y};

        /// A|| equilibrium value, used in corrector step
        fftw::mdbuffer<2u> aParEq_K{X, Y};

        /// Φ: the electrostatic potential.
        fftw::mdbuffer<2u> phi_K{X, Y}, phi_K_New{X, Y};

        /// sq(∇⊥) A∥, also parallel electron velocity
        fftw::mdbuffer<2u> ueKPar_K{X, Y}, ueKPar_K_New{X, Y};
        /// @}

        void for_each_xy(std::invocable<Dim, Dim> auto fun) {
            for (Dim x = 0; x < X; ++x) {
                for (Dim y = 0; y < Y; ++y) {
                    fun(x, y);
                }
            }
        }

        /// Iterate in phase space, will later be changed to account for phase space dims
        void for_each_kxky(std::invocable<Dim, Dim> auto fun) {
            for (Dim kx = 0; kx < X; ++kx) {
                for (Dim ky = 0; ky < Y; ++ky) {
                    fun(kx, ky);
                }
            }
        }

        void for_each_mxy(std::invocable<Dim, Dim, Dim> auto fun) {
            for (Dim m = 0; m < M; ++m) {
                for (Dim x = 0; x < X; ++x) {
                    for (Dim y = 0; y < Y; ++y) {
                        fun(m, x, y);
                    }
                }
            }
        }

        /// Returns a 2D mdspan of values in the XY space for a specified m.
        static auto sliceXY(fftw::mdbuffer<3u> &moments, Dim m) {
            return stdex::submdspan(moments.to_mdspan(), m, stdex::full_extent, stdex::full_extent);
        }

        /// Returns a DxDy of 2D mdspans for a specified m.
        static auto sliceXY(DxDy<fftw::mdbuffer<3u>> &moments, Dim m) {
            return DxDy{sliceXY(moments.DX, m), sliceXY(moments.DY, m)};
        }

        /// Prepares the δx and δy of viewPH in phase space
        void prepareDXY_PH(ViewXY viewPH, ViewXY viewPH_X, ViewXY viewPH_Y) {
            for_each_xy([&](Dim kx, Dim ky) {
                viewPH_X(kx, ky) = -kx_(kx) * 1i * viewPH(kx, ky);
                viewPH_Y(kx, ky) = -ky_(ky) * 1i * viewPH(kx, ky);
            });
        }


        /// computes bracket [view, other], also normalizing in the process
        void bracket(const ViewXY &dxOp1, const ViewXY &dyOp1, const ViewXY &dxOp2, const ViewXY &dyOp2,
                     const ViewXY &output) {
            for_each_xy([&](Dim x, Dim y) {
                output(x, y) = dxOp1(x, y) * dyOp2(x, y) - dyOp1(x, y) * dxOp2(x, y);
                output(x, y) /= double(X) * double(Y);
            });
        }

        /// bracket overload for DxDy params
        void bracket(const DxDy<ViewXY> &op1, const DxDy<ViewXY> &op2, const ViewXY &output) {
            bracket(op1.DX, op1.DY, op2.DX, op2.DY, output);
        }

        /// Bracket that only takes inputs and allocates temporaries and output
        [[nodiscard]] fftw::mdbuffer<2u> fullBracket(ViewXY op1, ViewXY op2);

        /// Compute derivatives in real space and store them in output
        void derivatives(const ViewXY &value, DxDy<ViewXY> output);

        /// Bracket that takes in derivatives that were already computed
        [[nodiscard]] fftw::mdbuffer<2u> halfBracket(DxDy<ViewXY> op1, DxDy<ViewXY> op2);

        // =================
        // Math helpers
        // TODO other file/class
        // =================

        [[nodiscard]] Real ky_(Dim ky) const { return Real(ky <= (KY / 2) ? ky : ky - KY) * Real(lx) / Real(ly); };
        [[nodiscard]] Real kx_(Dim kx) const { return Real(kx <= (KX / 2) ? kx : kx - KX); }; // TODO r2c

        [[nodiscard]] Real kPerp2(Dim kx, Dim ky) const {
            auto dkx = kx_(kx), dky = ky_(ky);
            return dkx * dkx + dky * dky;
        }

        [[nodiscard]] Real kPerp(Dim kx, Dim ky) const {
            return std::sqrt(kPerp2(kx, ky));
        }

        [[nodiscard]] Real exp_nu(Dim kx, Dim ky, Real niu2, Real dt) const {
            return std::exp(-(nu * kPerp2(kx, ky) + niu2 * std::pow(kPerp2(kx, ky), hyper_order)) * dt);
        }

        [[nodiscard]] Real exp_gm(Dim m, Real hyper_nuei, Real dt) const {
            return exp(-(Real(m) * nu_ei + std::pow(m, 2 * hyper_morder) * hyper_nuei) * dt);
        }

        [[nodiscard]] Real exp_eta(Dim kx, Dim ky, Real res2, Real dt) const {
            return std::exp(-(res * kPerp2(kx, ky) + res2 * std::pow(kPerp2(kx, ky), hyper_order)) * dt /
                            (1.0 + kPerp2(kx, ky) * de * de));
        }

        /// getTimestep calculates flows and magnetic fields to determine a dt.
        /// It also updates bPerpMax in the process.
        [[nodiscard]] Real getTimestep(DxDy<ViewXY> dPhi, DxDy<ViewXY> dNE, DxDy<ViewXY> dAPar) {
            // compute flows
            DxDy<ViewXY> ve, b;
            Real vyMax{0}, vxMax{0}, bxMax{0}, byMax{0};
            bPerpMax = 0;

            // Note that this is minus in Viriato, but we don't care because we're taking the absolute value anyway.
            ve.DX = dPhi.DY;
            ve.DY = dPhi.DX;
            b.DX = dAPar.DY;
            b.DY = dAPar.DX;

            for_each_xy([&](Dim x, Dim y) {
                bxMax = std::max(bxMax, std::abs(b.DX(x, y)));
                byMax = std::max(byMax, std::abs(b.DY(x, y)));
                bPerpMax = std::max(bPerpMax, std::sqrt(b.DX(x, y).real() * b.DX(x, y).real() +
                                                        b.DY(x, y).real() * b.DY(x, y).real()));
                vxMax = std::max(vxMax, std::abs(ve.DX(x, y)));
                vyMax = std::max(vyMax, std::abs(ve.DY(x, y)));
                if (rhoI >= smallRhoI) {
                    vxMax = std::max(vxMax, rhoS * rhoS * std::abs(dNE.DX(x, y)));
                    vyMax = std::max(vyMax, rhoS * rhoS * std::abs(dNE.DY(x, y)));
                }
            });

            Real kperpDum2 = std::pow(ky_(KY / 2), 2) + std::pow(Real(kx_(KX / 2)), 2);
            Real omegaKaw;
            if (rhoI < smallRhoI)
                omegaKaw = std::sqrt(1.0 + kperpDum2 * (3.0 / 4.0 * rhoI * rhoI + rhoS * rhoS))
                           * ky_(KY / 2 + 1) * bPerpMax / (1.0 + kperpDum2 * de * de);
            else
                omegaKaw = std::sqrt(kperpDum2 * (rhoS*rhoS - rhoI*rhoI / (Gamma0(0.5*kperpDum2*rhoI*rhoI)-1.0))) *
                    ky_(KY/2+1) * bPerpMax / std::sqrt(1.0 + kperpDum2 * de * de);

            Real dx = lx / Real(KX), dy = ly / Real(KY);
            Real CFLFlow = std::min({dx / vxMax, dy / vyMax, 2.0 / omegaKaw,
                                     std::min(dx / bxMax, dy / byMax) / rhoS / de / std::sqrt(M - 1)});

            // DEBUG
            std::cout << "CFLFlow: " << CFLFlow << std::endl;

            return CFLFrac * CFLFlow;
        }

        [[maybe_unused]] void print(ViewXY view) const {
            for (int x = 0; x < view.extent(0); ++x) {
                for (int y = 0; y < view.extent(1); ++y) {
                    std::cout << view(x, y) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "===========================================" << std::endl;
        }

        Real updateTimestep(Real dt, Real tempDt, bool noInc, Real relative_error) const;
    };
};