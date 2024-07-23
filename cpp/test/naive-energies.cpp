#include "util.hpp"
#include "Naive.h"
#include <gtest/gtest.h>

TEST(NaiveEnergy, OT01) {
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

TEST(NaiveEnergy, OT01NoDiffusion) {
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
    EXPECT_THAT(kin_final, AllClose(kin_init, 1e-5));
    EXPECT_THAT(mag_final + kin_final, AllClose(mag_init + kin_init, 1e-5));
}

TEST(NaiveEnergy, Gauss) {
    ahr::Naive naive{std::cout, 10, 16, 16};
    naive.init("gauss");
    auto [mag_init, kin_init] = naive.calculateEnergies();
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(kin_init, AllClose(0.0, 1e-7, 1e-6));

    naive.run(10, 0); // no saving
    auto [mag_final, kin_final] = naive.calculateEnergies();
    EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7));
    EXPECT_THAT(kin_final, AllClose(0.0, 1e-7, 1e-6));
    EXPECT_THAT(mag_final + kin_final, LeTolerant(mag_init + kin_init, 1e-7));

    // TODO find a good lower bound for energies here
}

TEST(NaiveEnergy, GaussNoDiffusion) {
    // turn off diffusion
    ahr::nu = 0;
    ahr::res = 0;

    ahr::Naive naive{std::cout, 10, 16, 16};
    naive.init("gauss");
    auto [mag_init, kin_init] = naive.calculateEnergies();

    naive.run(10, 0); // no saving
    auto [mag_final, kin_final] = naive.calculateEnergies();
    EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7));
    EXPECT_THAT(mag_final + kin_final, LeTolerant(mag_init + kin_init, 1e-7));

    EXPECT_THAT(mag_final, AllClose(mag_init, 1e-5));
    EXPECT_THAT(kin_final, AllClose(0.0, 1e-7, 1e-6));
    EXPECT_THAT(mag_final + kin_final, AllClose(mag_init + kin_init, 1e-5));
}

TEST(NaiveEnergy, GaussMagDiffusion) {
    // turn off kinetic diffusion
    ahr::nu = 0;
    ahr::res = 1.0;

    ahr::Naive naive{std::cout, 10, 16, 16};
    naive.init("gauss");
    auto [mag_init, kin_init] = naive.calculateEnergies();
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(kin_init, AllClose(0.0, 1e-7, 1e-6));

    naive.run(10, 0); // no saving
    auto [mag_final, kin_final] = naive.calculateEnergies();
    EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7));
    EXPECT_THAT(kin_final, AllClose(0.0, 1e-7, 1e-6));
    EXPECT_THAT(mag_final + kin_final, LeTolerant(mag_init + kin_init, 1e-7));

    // TODO good mag energy lower bound
}
