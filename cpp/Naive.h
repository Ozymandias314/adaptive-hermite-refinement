#pragma once

#include "HermiteRunner.h"

#include <fftw-cpp/fftw-cpp.h>

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
        fftw::plan<2u> plan{}, planInv{};
        Real initialDT{};

        /// \defgroup Buffers for all the physical quantities used.
        /// Names ending in PH mean the values are in phase space.
        /// @{

        /// g_m: moment values for moments $m \in [0,M-1]$.
        /// The following values are special moments:
        /// - m=0: n_e (charge density)
        /// - m=1: A|| (or Apar, parallel velocity)
        ///
        /// Additionally, moments is used only at the start and end.
        /// During the simulation, we use moments
        /// TODO maybe instead of these enormous amounts of memory, we could reuse (parallelism might suffer)
        fftw::mdbuffer<3u> moments, momentsPH, momentsPH_X, momentsPH_Y, moments_X, moments_Y;

        /// \phi: the electrostatic potential.
        fftw::mdbuffer<2u> phiPH;

        /// Temporary buffers used for various things.
        std::array<fftw::mdbuffer<2u>, 3> temp;

        /// @}

        static constexpr Dim A_PAR = 1;

        void for_each_xy(std::invocable<Dim, Dim> auto fun) {
            for (Dim x = 0; x < X; ++x) {
                for (Dim y = 0; y < Y; ++y) {
                    fun(x, y);
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

        /// sliceXY returns a 2D mdspan of values in the XY space for a specified m.
        static auto sliceXY(fftw::mdbuffer<3u> &moments, int m) {
            return stdex::submdspan(moments.to_mdspan(), m, stdex::full_extent, stdex::full_extent);
        }


        /// Computes bracket [viewPH, otherPH] assuming other has already been transformed into otherX and otherY.
        /// Uses other views for temporary storage, and stores results in viewPH_X.
        /// Currently really greedy.
        void bracketHalf(const ViewXY &viewPH, ViewXY viewPH_X, ViewXY viewPH_Y, ViewXY view_X, ViewXY view_Y,
                         const ViewXY &otherX,
                         const ViewXY &otherY, ViewXY bracket, ViewXY bracketPH) {
            for_each_xy([&](Dim kx, Dim ky) {
                viewPH_X(kx, ky) = -double(kx) * 1i * viewPH(kx, ky);
                viewPH_Y(kx, ky) = -double(ky) * 1i * viewPH(kx, ky);
            });

            planInv(viewPH_X, view_X);
            planInv(viewPH_Y, view_Y);

            for_each_xy([&](Dim x, Dim y) {
                bracket(x, y) = view_X(x, y) * otherY(x, y) - view_Y(x, y) * otherX(x, y);
            });

            plan(bracket, bracketPH);
        }
    };

}