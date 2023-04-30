#pragma once

#include "HermiteRunner.h"

#include <fftw-cpp/fftw-cpp.h>

namespace ahr {
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
        /// TODO reuse moments buffer for something else.
        fftw::mdbuffer<3u> moments, momentsPH;

        /// \phi: the electrostatic potential.
        fftw::mdbuffer<2u> phiPH;

        /// Temporary buffers used for various things.
        std::array<fftw::mdbuffer<2u>, 3> temp;

        /// @}

        void for_each_xy(std::invocable<Dim, Dim> auto fun) {
            for (Dim x = 0; x < X; ++x) {
                for (Dim y = 0; y < Y; ++y) {
                    fun(x, y);
                }
            }
        }

        void for_each_mxy(std::invocable<Dim, Dim, Dim> auto fun) {
            for (Dim m = 0; m < X; ++m) {
                for (Dim x = 0; x < X; ++x) {
                    for (Dim y = 0; y < Y; ++y) {
                        fun(m, x, y);
                    }
                }
            }
        }

        void mxy_copy(std::invocable<Dim, Dim,Dim > auto const &src, std::invocable<Dim, Dim,Dim > auto &dest) {
            for_each_mxy([&](Dim m, Dim x, Dim y) {
                dest(m, x, y) = src(m, x, y);
            });
        }

        /// sliceXY returns a 2D mdspan of values in the XY space for a specified m.
        static auto sliceXY(fftw::mdbuffer<3u> &moments, int m);
    };

}