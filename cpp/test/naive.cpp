#include "util.hpp"
#include "Naive.h"
#include <gtest/gtest.h>

TEST(BasicNaive, OT01Energy) {
    ahr::Naive naive{std::cout, 10, 16, 16};
    naive.init("OT01");
    auto [mag_init, kin_init] = naive.calculateEnergies();

    naive.run(10, 0); // no saving
    auto [mag_final, kin_final] = naive.calculateEnergies();
    EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7));
    EXPECT_THAT(kin_final, LeTolerant(kin_init, 1e-7));
    EXPECT_THAT(mag_final + kin_final, LeTolerant(mag_init + kin_init, 1e-7));

    // TODO find a good lower bound for energies here
}

TEST(BasicNaive, OT01NoDiffusionEnergy) {
    // turn off diffusion
    ahr::nu = 0;
    ahr::res = 0;

    ahr::Naive naive{std::cout, 10, 16, 16};
    naive.init("OT01");
    auto [mag_init, kin_init] = naive.calculateEnergies();

    naive.run(10, 0); // no saving
    auto [mag_final, kin_final] = naive.calculateEnergies();
    EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7));
    EXPECT_THAT(kin_final, LeTolerant(kin_init, 1e-7));
    EXPECT_THAT(mag_final + kin_final, LeTolerant(mag_init + kin_init, 1e-7));

    EXPECT_THAT(mag_final, AllClose(mag_init, 1e-5));
    EXPECT_THAT(kin_final, AllClose(mag_final, 1e-5));
    EXPECT_THAT(mag_final + kin_final, AllClose(mag_init + kin_init, 1e-5));
}