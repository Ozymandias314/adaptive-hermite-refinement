#include "Naive.h"

namespace ahr {
    namespace stdex = std::experimental;

    auto Naive::getMomentSlice(fftw::mdbuffer<3u> &moments, int m) {
        return stdex::submdspan(moments.to_mdspan(), m, stdex::full_extent, stdex::full_extent);
    }

    void Naive::init(mdspan<Real, dextents<Dim, 3u>> initialMoments, Dim N, Real initialDT) {
        this->N = N,
        M = initialMoments.extent(0),
        X = initialMoments.extent(1),
        Y = initialMoments.extent(2);
        assert(X == Y);

        fftw::mdbuffer<3u> moments{M, X, Y};
        for (int m = 0; m < M; ++m) {
            for (int x = 0; x < X; ++x) {
                for (int y = 0; y < Y; ++y) {
                    moments(m, x, y) = initialMoments(m, x, y);
                }
            }
        }

        auto inSubspan = getMomentSlice(moments, 0);
        auto outSubspan = getMomentSlice(momentsPH, 0);

        plan = fftw::plan<2u>::dft(inSubspan, outSubspan, fftw::FORWARD, fftw::MEASURE);

        for (int m = 0; m < M; ++m) {
            plan(getMomentSlice(moments, m), getMomentSlice(momentsPH, m));
        }
    }

    void Naive::run() {

        for (int t = 0; t < N; ++t) {
            // predictor step
            for (int m = 0; m < M; ++m) {
// TODO
//                plan(getMomentSlice(m), getMomentSlice());
            }

            // corrector step
            for (int m = 0; m < M; ++m) {

            }
        }
    }

    mdarray<Real, dextents<Dim, 3u>> Naive::getFinalValues() {
        mdarray<Real, dextents<Dim, 3u>> result{M, X, Y};
        for (int m = 0; m < M; ++m) {
            for (int x = 0; x < X; ++x) {
                for (int y = 0; y < Y; ++y) {
                    result(m, x, y) = momentsPH(m, x, y).real();
                }
            }
        }
        return result;
    }

    Naive::Naive(std::ostream &out) : HermiteRunner(out) {}
}