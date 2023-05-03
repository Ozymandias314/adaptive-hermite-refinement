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

            // Nabla
            for_each_xy([&](Dim kx, Dim ky) {
                auto dkx = double(kx), dky = double(ky), kPerp = dkx * dkx + dky * dky;
                // TODO minus or not? (kperp should already have a minus)
                nablaPerpAPar.PH(kx, ky) = kPerp * moments.PH(A_PAR, kx, ky);
            });

            // Compute N
            fullBracket(phi, sliceXY(moments, N_E));
            fullBracket(sliceXY(moments, A_PAR), nablaPerpAPar);
            // TODO use N

            // Compute A
            fullBracket(phi, sliceXY(moments, A_PAR));
            for_each_xy([&](Dim kx, Dim ky) {
                double const de = 1.2; //TODO
                nablaPerpAPar.PH(kx, ky) = de*de*nablaPerpAPar.PH(kx, ky);
            });

            // already using phi for temp storage
            fullBracket(phi, nablaPerpAPar, nablaPerpAPar.DX, nablaPerpAPar.PH_DX);
            // TODO use A

            // TODO compute G2

            for (int m = 3; m < M; ++m) {

            }

            // corrector step
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