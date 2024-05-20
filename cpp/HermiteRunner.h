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
         * @param N the number of timesteps.
         */
        virtual void
        init(Dim N) = 0;

        /**
         * run() will simulate Hermite moments for N timesteps.
         * init() must have been called before this call.
         */
        virtual void run(Dim saveInterval) = 0;

        /**
         *
         * @return Final values of APar.
         */
        virtual mdarray<Real, dextents<Dim, 2u>> getFinalAPar() = 0;


    protected:
        std::ostream &out;

    };

}