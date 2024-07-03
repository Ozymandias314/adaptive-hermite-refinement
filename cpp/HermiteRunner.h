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
         */
        virtual void
        init(std::string_view equilibriumName) = 0;

        /**
         * run() will simulate Hermite moments for N timesteps.
         * init() must have been called before this call.
         * @param N the number of timesteps.
         */
        virtual void run(Dim N, Dim saveInterval) = 0;

        /**
         *
         * @return Final values of APar.
         */
        virtual mdarray<Real, dextents<Dim, 2u>> getFinalAPar() = 0;


    protected:
        std::ostream &out;

    };

}