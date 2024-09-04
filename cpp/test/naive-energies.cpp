#include "Naive.h"
#include "NaiveTester.hpp"
#include "util.hpp"
#include <gtest/gtest.h>
#include <sstream>

using ahr::Real;

// Make this into a class if utilities need to be added
using NaiveEnergy = NaiveTester;

TEST_P(NaiveEnergy, Diffusion) {
  ahr::nu = 0.1;
  ahr::res = 0.1;
  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto [mag_init, kin_init] = naive.calculateEnergies();

  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(kin_init, AllClose(0.0, 1e-7, 1e-6));
  }

  naive.run(p.N, 0); // no saving
  auto [mag_final, kin_final] = naive.calculateEnergies();
  EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7, 1e-5));
  EXPECT_THAT(kin_final, LeTolerant(kin_init, 1e-7, 1e-5));
  EXPECT_THAT(mag_final + kin_final,
              LeTolerant(mag_init + kin_init, 1e-7, 1e-5));

  // TODO find a good lower bound for energies here
}

TEST_P(NaiveEnergy, NoDiffusion) {
  // turn off diffusion
  ahr::nu = 0;
  ahr::res = 0;

  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto [mag_init, kin_init] = naive.calculateEnergies();

  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(kin_init, AllClose(0.0, 1e-7, 1e-6));
  }

  naive.run(p.N, 0); // no saving
  auto [mag_final, kin_final] = naive.calculateEnergies();
  EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7, 1e-5));
  EXPECT_THAT(kin_final, LeTolerant(kin_init, 1e-7, 1e-5));
  EXPECT_THAT(mag_final + kin_final,
              LeTolerant(mag_init + kin_init, 1e-7, 1e-5));

  // The energy numerical error is roughly proportional to the # of timesteps
  auto rtol = 1e-6 * Real(p.N);
  EXPECT_THAT(mag_final, AllClose(mag_init, rtol, 1e-5));
  EXPECT_THAT(kin_final, AllClose(kin_init, rtol, 1e-5));
  EXPECT_THAT(mag_final + kin_final, AllClose(mag_init + kin_init, rtol, 1e-5));
}

TEST_P(NaiveEnergy, MagDiffusion) {
  // turn off kinetic diffusion
  ahr::nu = 0;
  ahr::res = 1.0;
  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto [mag_init, kin_init] = naive.calculateEnergies();
  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(kin_init, AllClose(0.0, 1e-7, 1e-5));
  }

  naive.run(p.N, 0); // no saving
  auto [mag_final, kin_final] = naive.calculateEnergies();

  // The energy numerical error is roughly proportional to the # of timesteps
  auto rtol = 1e-6 * Real(p.N);
  EXPECT_THAT(mag_final, LeTolerant(mag_init, 1e-7, 1e-5));
  EXPECT_THAT(kin_final, AllClose(kin_init, rtol, 1e-5));
  EXPECT_THAT(mag_final + kin_final,
              LeTolerant(mag_init + kin_init, 1e-7, 1e-5));
}

TEST_P(NaiveEnergy, KinDiffusion) {
  // turn off magnetic diffusion
  ahr::nu = 1.0;
  ahr::res = 0;
  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto [mag_init, kin_init] = naive.calculateEnergies();
  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(kin_init, AllClose(0.0, 1e-7, 1e-5));
  }

  naive.run(p.N, 0); // no saving
  auto [mag_final, kin_final] = naive.calculateEnergies();

  // The energy numerical error is roughly proportional to the # of timesteps
  auto rtol = 1e-6 * Real(p.N);
  EXPECT_THAT(mag_final, AllClose(mag_init, rtol, 1e-5));
  EXPECT_THAT(kin_final, LeTolerant(kin_init, 1e-7, 1e-5));
  EXPECT_THAT(mag_final + kin_final,
              LeTolerant(mag_init + kin_init, 1e-7, 1e-5));
}

using namespace testing;
INSTANTIATE_TEST_SUITE_P(
    NaiveEnergyTestsSmallM, NaiveEnergy,
    ConvertGenerator<TesterParam::Tuple>(Combine(Values(2, 4),       // M
                                                 Values(16, 32, 64), // X
                                                 Values(10, 20),     // N
                                                 Values("gauss",
                                                        "OT01")) // eq
                                         ));

INSTANTIATE_TEST_SUITE_P(
    NaiveEnergyTestsLargeM, NaiveEnergy,
    ConvertGenerator<TesterParam::Tuple>(Combine(Values(10, 20),          // M
                                                 Values(16),              // X
                                                 Values(10, 20),          // N
                                                 Values("gauss", "OT01")) // eq
                                         ));