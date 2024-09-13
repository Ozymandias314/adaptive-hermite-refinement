#pragma once
#include "Naive.h"
#include <gtest/gtest.h>

struct TesterParam {
  ahr::Dim M, X, N;
  std::string_view eq; ///< equillibrium

  /// M, X, N, eq
  using Tuple = std::tuple<ahr::Dim, ahr::Dim, ahr::Dim, std::string_view>;

  explicit TesterParam(Tuple t)
      : M(std::get<0>(t)), X(std::get<1>(t)), N(std::get<2>(t)),
        eq(std::get<3>(t)) {}
};

struct TesterParam2 {
  ahr::Dim M, X, N;
  ahr::Real nu, res;

  /// M, X, N, nu, res
  using Tuple = std::tuple<ahr::Dim, ahr::Dim, ahr::Dim, ahr::Real, ahr::Real>;

  explicit TesterParam2(Tuple t)
      : M(std::get<0>(t)), X(std::get<1>(t)), N(std::get<2>(t)),
        nu(std::get<3>(t)), res(std::get<4>(t)) {}
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

class NaiveTester : public TesterWithOutput,
                    public ::testing::WithParamInterface<TesterParam::Tuple> {};
class NaiveTester2 : public TesterWithOutput,
                     public ::testing::WithParamInterface<TesterParam2::Tuple> {
};