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
        fft = fftw::plan<2u>::dft(sliceXY(momentsReal, 0), sliceXY(moments.PH, 0), fftw::FORWARD, fftw::MEASURE);
        fftInv = fftw::plan<2u>::dft(sliceXY(moments.PH, 0), sliceXY(momentsReal, 0), fftw::BACKWARD, fftw::MEASURE);

        // Initialize moments
        for_each_mxy([&](Dim m, Dim x, Dim y) {
            momentsReal(m, x, y) = initialMoments(m, x, y);
        });

        // Transform moments into phase space
        for (int m = 0; m < M; ++m) {
            fft(sliceXY(momentsReal, m), sliceXY(moments.PH, m));
        }
    }

    void Naive::run() {

        for (int t = 0; t < N; ++t) {
            // predictor step

            // Phi
            // TODO maybe compute PHI? Or maybe that happened in the previous timestep
            prepareDXY_PH(phi.PH, phi.PH_DX, phi.PH_DY);
            fftInv(phi.PH_DX, phi.DX);
            fftInv(phi.PH_DY, phi.DY);

            // Nabla
            for_each_xy([&](Dim kx, Dim ky) {
                auto dkx = double(kx), dky = double(ky), kPerp = dkx * dkx + dky * dky;
                // TODO minus or not? (kperp should already have a minus)
                nablaPerpAPar.PH_DX(kx, ky) = -dkx * 1i * kPerp * nablaPerpAPar.PH(kx, ky);
                nablaPerpAPar.PH_DY(kx, ky) = -dky * 1i * kPerp * nablaPerpAPar.PH(kx, ky);
            });
            fftInv(nablaPerpAPar.PH_DX, nablaPerpAPar.DX);
            fftInv(nablaPerpAPar.PH_DY, nablaPerpAPar.DY);

            for (int m = 2; m < M; ++m) {
                ViewXY viewPH_X = sliceXY(moments.PH_DX, m), viewPH_Y = sliceXY(moments.PH_DY, m);
                ViewXY viewX = sliceXY(moments.DX, m), viewY = sliceXY(moments.DY, m);

                prepareDXY_PH(sliceXY(moments.PH, m), viewPH_X, viewPH.DY);
                fftInv(viewPH_X, viewX);
                fftInv(viewPH_Y, viewY);

                ViewXY otherX = sliceXY(moments.DX, A_PAR),
                        otherY = sliceXY(moments.DY, A_PAR);
                bracket(viewX, viewY, otherX, otherY, temp[0]);

                fft(temp[0], temp[1]);
            }

            // corrector step
            for (int m = 0; m < M; ++m) {

            }
        }
    }

    mdarray<Real, dextents<Dim, 3u>> Naive::getFinalValues() {
        for (int m = 0; m < M; ++m) {
            fftInv(sliceXY(moments.PH, m), sliceXY(momentsReal, m));
        }

        mdarray<Real, dextents<Dim, 3u>> result{M, X, Y};
        for_each_mxy([&](Dim m, Dim x, Dim y) {
            result(m, x, y) = momentsReal(m, x, y).real();
        });

        return result;
    }
}