#include "Naive.h"
#include "NaiveTester.hpp"
#include "util.hpp"
#include <gtest/gtest.h>
#include <sstream>

using ahr::Real;

// Make this into a class if utilities need to be added
using NaiveEnergy = NaiveTester;
using Energies = ahr::Naive::Energies;

Energies expectedEnergies(Real t, Energies e_init) {
  return {.magnetic = std::exp(-t * ahr::res * 2) * e_init.magnetic,
          .kinetic = std::exp(-t * ahr::nu * 2) * e_init.kinetic};
}

#define CHECK_ENERGIES()

TEST_P(NaiveEnergy, Diffusion) {
  ahr::nu = 0.1;
  ahr::res = 0.1;
  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto const e_init = naive.calculateEnergies();

  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(e_init.kinetic, AllClose(0.0, 1e-7, 1e-6));
  }

  naive.run(p.N, 0); // no saving
  auto const e_final = naive.calculateEnergies();
  EXPECT_THAT(e_final.magnetic, LeTolerant(e_init.magnetic, 1e-7, 1e-5));
  EXPECT_THAT(e_final.kinetic, LeTolerant(e_init.kinetic, 1e-7, 1e-5));
  EXPECT_THAT(e_final.magnetic + e_final.kinetic,
              LeTolerant(e_init.magnetic + e_init.kinetic, 1e-7, 1e-5));

  // TODO find a good lower bound for energies here
}

TEST_P(NaiveEnergy, NoDiffusion) {
  // turn off diffusion
  ahr::nu = 0;
  ahr::res = 0;

  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto const e_init = naive.calculateEnergies();

  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(e_init.kinetic, AllClose(0.0, 1e-7, 1e-6));
  }

  naive.run(p.N, 0); // no saving
  auto const e_final = naive.calculateEnergies();
  EXPECT_THAT(e_final.magnetic, LeTolerant(e_init.magnetic, 1e-7, 1e-5));
  EXPECT_THAT(e_final.kinetic, LeTolerant(e_init.kinetic, 1e-7, 1e-5));
  EXPECT_THAT(e_final.magnetic + e_final.kinetic,
              LeTolerant(e_init.magnetic + e_init.kinetic, 1e-7, 1e-5));

  // The energy numerical error is roughly proportional to the # of timesteps
  auto rtol = 1e-6 * Real(p.N);
  EXPECT_THAT(e_final.magnetic, AllClose(e_init.magnetic, rtol, 1e-5));
  EXPECT_THAT(e_final.kinetic, AllClose(e_init.kinetic, rtol, 1e-5));
  EXPECT_THAT(e_final.magnetic + e_final.kinetic,
              AllClose(e_init.magnetic + e_init.kinetic, rtol, 1e-5));
}

TEST_P(NaiveEnergy, MagDiffusion) {
  // turn off kinetic diffusion
  ahr::nu = 0;
  ahr::res = 1.0;
  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto const e_init = naive.calculateEnergies();
  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(e_init.kinetic, AllClose(0.0, 1e-7, 1e-5));
  }

  naive.run(p.N, 0); // no saving
  auto const e_final = naive.calculateEnergies();

  // The energy numerical error is roughly proportional to the # of timesteps
  auto const rtol = 1e-6 * Real(p.N);
  EXPECT_THAT(e_final.magnetic, LeTolerant(e_init.magnetic, 1e-7, 1e-5));
  EXPECT_THAT(e_final.kinetic, AllClose(e_init.kinetic, rtol, 1e-5));
  EXPECT_THAT(e_final.magnetic + e_final.kinetic,
              LeTolerant(e_init.magnetic + e_init.kinetic, 1e-7, 1e-5));
}

TEST_P(NaiveEnergy, KinDiffusion) {
  // turn off magnetic diffusion
  ahr::nu = 1.0;
  ahr::res = 0;
  auto p = TesterParam{GetParam()};

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init(p.eq);
  auto const e_init = naive.calculateEnergies();
  if (p.eq == "gauss") {
    // Kinetic energy should be zero (absolute tolerance hence needed)
    EXPECT_THAT(e_init.kinetic, AllClose(0.0, 1e-7, 1e-5));
  }

  naive.run(p.N, 0); // no saving
  auto e_final = naive.calculateEnergies();

  // The energy numerical error is roughly proportional to the # of timesteps
  auto rtol = 1e-6 * Real(p.N);
  EXPECT_THAT(e_final.magnetic, AllClose(e_init.magnetic, rtol, 1e-5));
  EXPECT_THAT(e_final.kinetic, LeTolerant(e_init.kinetic, 1e-7, 1e-5));
  EXPECT_THAT(e_final.magnetic + e_final.kinetic,
              LeTolerant(e_init.magnetic + e_init.kinetic, 1e-7, 1e-5));
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