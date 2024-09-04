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

class NaiveTester : public ::testing::TestWithParam<TesterParam::Tuple> {
protected:
  // Ignore output
  std::ostringstream out;

  void TearDown() override {
    if (HasFailure()) {
      std::cout << "Full output:" << out.str() << std::endl;
    }
  }
};