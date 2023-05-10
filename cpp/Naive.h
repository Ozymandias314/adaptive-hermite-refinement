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

        void init(mdspan<Real, dextents<Dim, 3u>> initialMoments, Dim N_, Real initialDT_) override;

        void run() override;

        mdarray<Real, dextents<Dim, 3u>> getFinalValues() override;

        using ViewXY = stdex::mdspan<Complex, stdex::dextents<Dim, 2u>>;
    private:
        Dim M, X, Y, N{};
        fftw::plan<2u> fft{}, fftInv{};
        Real initialDT{};

        static constexpr Dim N_E = 0;
        static constexpr Dim A_PAR = 1;

        static constexpr Dim N_TEMP_BUFFERS = 3;

        template<class Buffer>
        struct BracketBuf {
            using buffer_t = Buffer;
            Buffer PH, PH_DX, PH_DY, DX, DY;

            template<typename ...Args>
            requires std::constructible_from<Buffer, Args...>
            explicit BracketBuf(Args &&...args) :
                    PH{std::forward<Args>(args)...},
                    PH_DX{std::forward<Args>(args)...},
                    PH_DY{std::forward<Args>(args)...},
                    DX{std::forward<Args>(args)...},
                    DY{std::forward<Args>(args)...} {}

            BracketBuf(Buffer ph, Buffer phDx, Buffer phDy, Buffer dx, Buffer dy)
                    : PH(ph), PH_DX(phDx), PH_DY(phDy), DX(dx), DY(dy) {}

            template<class U>
            requires std::convertible_to<Buffer, U>
            operator BracketBuf<U>() { // NOLINT(google-explicit-constructor)
                return {
                        U(PH), U(PH_DX), U(PH_DY), U(DX), U(DY)
                };
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
        fftw::mdbuffer<3u> moments{M, X, Y};
        fftw::mdbuffer<3u> moments_K{M, X, Y};

        /// Φ: the electrostatic potential.
        fftw::mdbuffer<2u> phi_K{X, Y};

        /// sq(∇⊥) A∥
        fftw::mdbuffer<2u> nablaPerpAPar_K{X, Y};

        /// Temporary buffers used for various things.
        std::array<fftw::mdbuffer<2u>, N_TEMP_BUFFERS> temp;

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
        static auto sliceXY(fftw::mdbuffer<3u> &moments, int m) {
            return stdex::submdspan(moments.to_mdspan(), m, stdex::full_extent, stdex::full_extent);
        }

        /// Returns a BracketBuf of 2D mdspans for a specified m.
        static auto sliceXY(BracketBuf<fftw::mdbuffer<3u>> &moments, int m) {
            return BracketBuf{
                    sliceXY(moments.PH, m),
                    sliceXY(moments.PH_DX, m),
                    sliceXY(moments.PH_DY, m),
                    sliceXY(moments.DX, m),
                    sliceXY(moments.DY, m)
            };
        }

        /// Prepares the δx and δy of viewPH in phase space
        void prepareDXY_PH(ViewXY viewPH, ViewXY viewPH_X, ViewXY viewPH_Y) {
            for_each_xy([&](Dim kx, Dim ky) {
                viewPH_X(kx, ky) = -double(kx) * 1i * viewPH(kx, ky);
                viewPH_Y(kx, ky) = -double(ky) * 1i * viewPH(kx, ky);
            });
        }


        /// computes bracket [view, other]
        void bracket(const ViewXY &viewX, const ViewXY &view_Y, const ViewXY &otherX, const ViewXY &otherY,
                     const ViewXY &bracket) {
            for_each_xy([&](Dim x, Dim y) {
                bracket(x, y) = viewX(x, y) * otherY(x, y) - view_Y(x, y) * otherX(x, y);
            });
        }

        void fullBracket(BracketBuf<ViewXY> op1, BracketBuf<ViewXY> op2, ViewXY brack, ViewXY brackPH) {
            prepareDXY_PH(op1.PH, op1.PH_DX, op1.PH_DY);
            fftInv(op1.PH_DX, op1.DX);
            fftInv(op1.PH_DY, op1.DY);
            prepareDXY_PH(op2.PH, op2.PH_DX, op2.PH_DY);
            fftInv(op2.PH_DX, op2.DX);
            fftInv(op2.PH_DY, op2.DY);

            bracket(op1.DX, op1.DY, op2.DX, op2.DY, brack);
            fft(brack, brackPH);
        }

        void fullBracket(BracketBuf<ViewXY> op1, BracketBuf<ViewXY> op2) {
            fullBracket(op1, op2, op1.DX, op1.PH_DX);
        };

        /// Bracket that only takes inputs and allocates temporaries and output
        [[nodiscard]] fftw::mdbuffer<2u> fullBracket(ViewXY op1, ViewXY op2);;

        // =================
        // Math helpers
        // TODO other file/class
        // =================

        // Nonlinear operators


        void NonlinearN(const ViewXY &bracketPhiNE_K, const ViewXY &bracketAParNablaPerpAPar_K, ViewXY result) {
            for_each_kxky([&](Dim kx, Dim ky) {
                result(kx, ky) = -bracketPhiNE_K(kx, ky) + bracketAParNablaPerpAPar_K(kx, ky);
            });
        }

        void NonlinearA(const ViewXY &bracketPhiAPar_K, const ViewXY &bracketPhiDeNablaPerpAPar_K,
                        const ViewXY &bracketNeG2APar_K, ViewXY result) {
            for_each_kxky([&](Dim kx, Dim ky) {
                result(kx, ky) = (-bracketPhiAPar_K(kx, ky) + bracketPhiDeNablaPerpAPar_K(kx, ky) +
                                  rhoS * rhoS * bracketNeG2APar_K(kx, ky)) / (1 + de * de * kPerp(kx, ky));
            });
        }

        void NonlinearG2(const ViewXY &bracketPhiG2_K, const ViewXY &bracketAParG3_K,
                         const ViewXY &bracketAParNablaPerpAPar_K, ViewXY result) {
            for_each_kxky([&](Dim kx, Dim ky) {
                result(kx, ky) = -bracketPhiG2_K(kx, ky) +
                                 std::sqrt(3.0) * rhoS / de * bracketAParG3_K(kx, ky) +
                                 std::sqrt(2.0) * bracketAParNablaPerpAPar_K(kx, ky);
            });
        }

        void NonlinearGM(Dim m, const ViewXY &bracketPhiGM_K, const ViewXY &bracketAParGMMinus_K,
                         const ViewXY &bracketAParGMPlus_K, ViewXY result) {
            for_each_kxky([&](Dim kx, Dim ky) {
                result(kx, ky) = -bracketPhiGM_K(kx, ky) +
                                 std::sqrt(m - 1) * rhoS / de * bracketAParGMMinus_K(kx, ky) +
                                 std::sqrt(m) * rhoS / de * bracketAParGMPlus_K(kx, ky);
            });
        }
    };
};