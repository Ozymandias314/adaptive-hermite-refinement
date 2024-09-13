#include "Naive.h"
#include "NaiveTester.hpp"
#include "util.hpp"
#include <gtest/gtest.h>
#include <sstream>

using ahr::Real;
using namespace ::testing;

// Make this into a class if utilities need to be added
using NaiveEnergy = NaiveTester;
using NaiveEnergy2 = NaiveTester2;
using Energies = ahr::Naive::Energies;

Energies expectedEnergies(Real t, Energies e_init) {
  return {.magnetic = std::exp(-t * ahr::res * 2) * e_init.magnetic,
          .kinetic = std::exp(-t * ahr::nu * 2) * e_init.kinetic};
}

#define CHECK_ENERGIES()

// TODO
// - gauss:
//    - start at 32, kinetic 0, magnetic decreases with diffusion
// - OT01:
//    - only look at sum
//    - check that sum diffusion is roughly the sum of the expected diffusions
//    - total energy shouldn't decrease
TEST_P(NaiveEnergy2, Gauss) {
  auto p = TesterParam2{GetParam()};
  ahr::nu = p.nu;
  ahr::res = p.res;

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init("gauss");
  auto const e_init = naive.calculateEnergies();

  // Kinetic energy should be zero (absolute tolerance hence needed)
  EXPECT_THAT(e_init.kinetic, AllClose(0.0, 1e-7, 1e-6));

  naive.run(p.N, 0); // no saving
  auto const e_final = naive.calculateEnergies();
  auto const e_expected = expectedEnergies(naive.elapsedTime(), e_init);

  // Kinetic should still be zero
  EXPECT_THAT(e_final.kinetic, AllClose(0.0, 1e-7, 1e-5));
  EXPECT_THAT(e_expected.kinetic, AllClose(0.0, 1e-7, 1e-5));

  // Magnetic energy should not increase
  EXPECT_THAT(e_final.magnetic, LeTolerant(e_init.magnetic, 1e-7, 1e-5));

  // Magnetic energy should roughly diffuse.
  Real diffused = e_init.magnetic - e_final.magnetic;
  Real expected_diffusion = e_init.magnetic - e_expected.magnetic;

  // Add energy tolerance
  auto const rtol = 2e-5 * Real(p.N);
  auto const e_tol = rtol * e_init.magnetic;
  EXPECT_THAT(diffused,
              AllOf(GeTolerant(expected_diffusion * 0.9, 0.0, e_tol),
                    LeTolerant(expected_diffusion * 1.1, 0.0, e_tol)));

  std::cout << "mag_init: " << e_init.magnetic
            << " mag_final: " << e_final.magnetic
            << " mag_expected: " << e_expected.magnetic << std::endl;
  std::cout << "diffused: " << diffused
            << " expected_diffusion: " << expected_diffusion << std::endl;
}

TEST_P(NaiveEnergy2, OT01) {
  auto p = TesterParam2{GetParam()};
  ahr::nu = p.nu;
  ahr::res = p.res;

  ahr::Naive naive{out, p.M, p.X, p.X};
  naive.init("OT01");
  auto const e_init = naive.calculateEnergies();

  naive.run(p.N, 0); // no saving
  auto const e_final = naive.calculateEnergies();
  auto const e_expected = expectedEnergies(naive.elapsedTime(), e_init);

  // Energy sum should not increase
  EXPECT_THAT(e_final.total(), LeTolerant(e_init.total(), 1e-7, 1e-5));

  // Energy sum should roughly diffuse.
  Real const diffused = e_init.total() - e_final.total();
  Real const expected_diffusion = e_init.total() - e_expected.total();

  // Add energy tolerance
  auto const rtol = 2e-5 * Real(p.N);
  auto const e_tol = rtol * e_init.total();
  EXPECT_THAT(diffused,
              AllOf(GeTolerant(expected_diffusion * 0.9, 0.0, e_tol),
                    LeTolerant(expected_diffusion * 1.1, 0.0, e_tol)));

  std::cout << "mag_init: " << e_init.magnetic
            << " mag_final: " << e_final.magnetic
            << " mag_expected: " << e_expected.magnetic << std::endl;
  std::cout << "kin_init: " << e_init.kinetic
            << " kin_final: " << e_final.kinetic
            << " kin_expected: " << e_expected.kinetic << std::endl;
  std::cout << "diffused: " << diffused
            << " expected_diffusion: " << expected_diffusion << std::endl;
}

using namespace testing;
INSTANTIATE_TEST_SUITE_P(
    NaiveEnergy2TestsSmallM, NaiveEnergy2,
    ConvertGenerator<TesterParam2::Tuple>(
        Combine(Values(2, 4),          // M
                Values(32, 64),        // X
                Values(10, 20),        // N
                Values(0.0, 0.1, 1.0), // nu
                Values(0.0) //, 0.1, 1.0) // res - TODO magnetic diffusion
                )));
INSTANTIATE_TEST_SUITE_P(
    NaiveEnergy2TestsLargeM, NaiveEnergy2,
    ConvertGenerator<TesterParam2::Tuple>(
        Combine(Values(10, 20),        // M
                Values(32),            // X
                Values(10, 20),        // N
                Values(0.0, 0.1, 1.0), // nu
                Values(0.0) //, 0.1, 1.0)) // res - TODO magnetic diffusion
                )));

// TODO remove old tests
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

  // Expected energy diffusion
  auto const e_expected = expectedEnergies(naive.elapsedTime(), e_init);
  auto const rtol = 1e-6 * Real(p.N);

  //  TODO magnetic diffusion (or the estimate) is broken
  //  EXPECT_THAT(e_final.magnetic, AllClose(e_expected.magnetic, rtol, 1e-5));
  EXPECT_THAT(e_final.kinetic, AllClose(e_expected.kinetic, rtol, 1e-5));
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

  // Estimated energy diffusion
  auto const e_expected = expectedEnergies(naive.elapsedTime(), e_init);
  //  TODO magnetic diffusion (or the estimate) is broken
  //  EXPECT_THAT(e_final.magnetic, AllClose(e_expected.magnetic, rtol, 1e-5));
  EXPECT_THAT(e_final.kinetic, AllClose(e_expected.kinetic, rtol, 1e-5));
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
  EXPECT_THAT(e_final.kinetic, LeTolerant(e_init.kinetic, 1e-7, 1e-5));

  // Expected energy diffusion
  auto const e_expected = expectedEnergies(naive.elapsedTime(), e_init);
  EXPECT_THAT(e_final.magnetic, AllClose(e_expected.magnetic, rtol, 1e-5));
  EXPECT_THAT(e_final.kinetic, AllClose(e_expected.kinetic, rtol, 1e-5));
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