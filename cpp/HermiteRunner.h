#pragma once

#include "typedefs.h"
#include <experimental/mdspan>
#include <experimental/mdarray>
#include <iostream>

namespace ahr {
    using std::experimental::mdspan;
    using std::experimental::mdarray;
    using std::experimental::dextents;

    class HermiteRunner {
    public:

        explicit HermiteRunner(std::ostream &out);

        /**
         * init prepares the hermite simulation.
         * @param initialMoments the initial values of the moments. The dimensions of this value also encode the
         * dimensions (M, X, Y, in this order) of the simulation.
         * @param N the number of timesteps.
         * @param initialDT The starting integration timestep.
         */
        virtual void init(mdspan<Real, dextents<Dim, 3u>> initialMoments, Dim N, Real initialDT) = 0;

        /**
         * run() will simulate Hermite moments for N timesteps.
         * init() must have been called before this call.
         */
        virtual void run() = 0;

        /**
         *
         * @return Final values of the moments.
         */
        virtual mdarray<Real, dextents<Dim, 3u>> getFinalValues() = 0;


    protected:
        std::ostream &out;

    };

}