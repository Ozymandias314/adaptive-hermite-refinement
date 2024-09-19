#pragma once
#include "Naive.h"
#include <gtest/gtest.h>

namespace ahr {

struct TesterParam {
  ahr::Dim M, X, N;
  ahr::Real nu, res;

  /// M, X, N, nu, res
  using Tuple = std::tuple<ahr::Dim, ahr::Dim, ahr::Dim, ahr::Real, ahr::Real>;

  explicit TesterParam(Tuple t)
      : M(std::get<0>(t)), X(std::get<1>(t)), N(std::get<2>(t)),
        nu(std::get<3>(t)), res(std::get<4>(t)) {}

  static std::string d_to_str(double d) {
    auto r = std::to_string(std::round(d * 100) / 100);
    std::ranges::replace(r, '.', '_');
    return r;
  }

  std::string nu_str() const { return d_to_str(nu); }

  std::string res_str() const { return d_to_str(res); }

  std::string to_param_str() const {
    std::ostringstream oss;
    oss << "M" << M << "_X" << X << "_N" << N << "_nu" << nu_str() << "_res"
        << res_str();

    return oss.str();
  }
};

class TesterWithOutput : public ::testing::Test {
protected:
  // Ignore output
  std::ostringstream out{};

  void TearDown() override {
    if (HasFailure()) {
      std::cout << "Full output:" << out.str() << std::endl;
    }
  }
};

template <class Param = TesterParam>
class NaiveTester : public TesterWithOutput,
                    public testing::WithParamInterface<Param> {
public:
  class Printer {
  public:
    std::string
    operator()(const testing::TestParamInfo<Param> &info) const {
      return info.param.to_param_str();
    }
  };
};

} // namespace ahr
