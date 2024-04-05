#pragma once

#include "constants.h"
#include "Naive.h"
#include <experimental/mdarray>
#include <utility>

namespace ahr {
    namespace stdex = std::experimental;

    // Note that these differ from julia because they're 0-indexed.
    Real xx(Dim x, Dim X) {
        return lx * (Real(x) - Real(X) / 2.0) / Real(X);
    }

    Real yy(Dim y, Dim Y) {
        return ly * (Real(y) - Real(Y) / 2.0) / Real(Y);
    }


    auto equilibriumGauss(Dim X, Dim Y) {
        Naive::Buf2D aParEq{X, Y}, phiEq{X, Y};

        for (int x = 0; x < X; ++x) {
            for (int y = 0; y < Y; ++y) {
                using std::numbers::pi;
                aParEq(x, y) = a0 * std::exp(-std::pow(yy(y, Y) * 2 * pi * 2 / ly, 2)) *
                               std::exp(-std::pow(xx(x, X) * 2 * pi * 2 / lx, 2));

                phiEq(x, y) = 0;
            }
        }
        return std::make_pair(std::move(aParEq), std::move(phiEq));
    }

    auto equilibriumOT01(Dim X, Dim Y) {
        Naive::Buf2D aParEq{X, Y}, phiEq{X, Y};

        for (int x = 0; x < X; ++x) {
            for (int y = 0; y < Y; ++y) {
                using std::numbers::pi;
                aParEq(x, y) = std::cos(4*pi*xx(x, X)/lx)+2*std::cos(2*pi*yy(y, Y)/ly);
                phiEq(x, y) = -2 * (std::cos(2 * pi * xx(x, X) / lx) + std::cos(2 * pi * yy(y, Y) / ly));
            }
        }
        return std::make_pair(std::move(aParEq), std::move(phiEq));
    }
}
