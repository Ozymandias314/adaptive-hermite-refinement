#include "Naive.h"
#include "NaiveTester.hpp"
#include "util.hpp"

#include <cnpy.h>

#include <filesystem>
#include <gtest/gtest.h>

namespace ahr {

// TODO WithEquilibrium
struct MomentParam : TesterParam {
  std::string equilibrium;

  /// eq, M, X, N, nu, res
  using Tuple = std::tuple<std::string, Dim, Dim, Dim, Real, Real>;

  explicit MomentParam(Tuple t)
      : equilibrium(std::get<0>(t)), TesterParam{slice<1>(t)} {}

  std::string to_param_str() const {
    return equilibrium + "_" + TesterParam::to_param_str();
  }
};

class NaiveMoments : public NaiveTester<MomentParam> {
protected:
  std::string getFilename(Dim m) {
    std::ostringstream oss;
    auto const p = GetParam();

    // Tests are run from inside the `test` directory
    oss << "./_test_data/" << p.equilibrium << "_M" << p.M << "_x" << p.X
        << "_n" << p.N << "_res" << p.res_str() << "_nu" << p.nu_str() << "_m"
        << m << ".npy";

    return oss.str();
  }

  // Owning holder of a npy array with convenience to see it as an mdspan.
  // We can use this to avoid needlessly copying into an mdarray.
  class NpyMdspan {
    cnpy::NpyArray array_;

  public:
    explicit NpyMdspan(const cnpy::NpyArray &array) : array_(array) {}

    // TODO(luka) const view
    Naive::ViewXY view() {
      std::span<size_t, 2> extents{array_.shape.data(), 2};
      return Naive::ViewXY{array_.data<Real>(), extents};
    }

    [[nodiscard]] bool valid() const {
      return array_.word_size == sizeof(Real);
    }
  };

  NpyMdspan readMoment(Dim m) {
    auto const filename = getFilename(m);
    return NpyMdspan{cnpy::npy_load(filename)};
  }
};

TEST_P(NaiveMoments, CheckMoments) {
  auto const p = GetParam();

  nu = p.nu;
  res = p.res;

  Naive naive{out, p.M, p.X, p.X};

  naive.init(p.equilibrium);
  naive.run(p.N + 1, 0); // no saving

  for (Dim m = 0; m < p.M; m++) {
    // To update values, uncomment these 2 lines
    // std::cout << "WARNING!: Overwriting " << getFilename(m) << std::endl;
    // naive.exportToNpy(getFilename(m), naive.getMoment(m));

    auto npy = readMoment(m);
    ASSERT_TRUE(npy.valid());

    EXPECT_THAT(naive.getMoment(m).to_mdspan(),
                MdspanElementsAllClose(npy.view(), 1e-14, 1e-15))
        << "Moment " << m << " mismatch!";
  }
}

using namespace testing;
INSTANTIATE_TEST_SUITE_P(NaiveMomentsSmallM, NaiveMoments,
                         ConvertGenerator<MomentParam::Tuple>(Combine(
                             Values("OT01", "gauss"), // equilibrium
                             Values(2, 4),            // M
                             Values(32, 64, 128),     // X
                             Values(20),              // N
                             Values(0.1, 1.0),     // nu - if 0, energy blows up
                             Values(0.0, 0.1, 1.0) // res
                             )),
                         NaiveMoments::Printer{});

INSTANTIATE_TEST_SUITE_P(NaiveMomentsLargeM, NaiveMoments,
                         ConvertGenerator<MomentParam::Tuple>(Combine(
                             Values("OT01", "gauss"), // equilibrium
                             Values(10, 20),          // M
                             Values(32, 64),          // X
                             Values(20),              // N
                             Values(0.1, 1.0), // nu - if 0, energy blows up
                             Values(0.1, 1.0)  // res
                             )),
                         NaiveMoments::Printer{});

}; // namespace ahr