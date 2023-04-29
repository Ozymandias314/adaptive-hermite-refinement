#pragma once

#include "HermiteRunner.h"

#include <fftw-cpp/fftw-cpp.h>

namespace ahr {
    class Naive : public ahr::HermiteRunner {
    public:
        explicit Naive(std::ostream &out);

        void init(mdspan<Real, dextents<Dim, 3u>> initialMoments, Dim N, Real initialDT) override;

        void run() override;

        mdarray<Real, dextents<Dim, 3u>> getFinalValues() override;

    private:
        fftw::mdbuffer<3u> momentsPH{0, 0, 0};
        Dim M{}, X{}, Y{}, N{};
        fftw::plan<2u> plan{};

        auto getMomentSlice(fftw::mdbuffer<3u> &moments, int m);
    };

}