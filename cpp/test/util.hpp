#pragma once

#include <fftw-cpp/fftw-cpp.h>
#include <gmock/gmock.h>

using ::testing::PrintToString;
MATCHER_P3(AllClose, val, rel_tol, abs_tol,
           PrintToString(val) + " ±" + PrintToString(abs_tol) + " (±" + PrintToString(double(rel_tol) * std::abs(val)) +
           ")") {
    auto diff = std::max(arg, val) - std::min(arg, val);
    double tolerance_diff = double(diff) - double(abs_tol) - double(rel_tol) * std::abs(val);


    return tolerance_diff <= 0;
}


// abs is zero by default
template<typename V, typename R>
auto AllClose(V &&val, R &&rel_tol) {
    return AllClose(std::forward<V>(val), std::forward<R>(rel_tol), 0.0);
}

MATCHER_P3(LeTolerant, val, rel_tol, abs_tol,
           "less or equal than " + PrintToString(val) + " ±" + PrintToString(abs_tol) + " (±" +
           PrintToString(double(rel_tol) * std::abs(val)) + ")") {
    auto diff = double(arg) - double(val);
    double tolerance_diff = double(diff) - double(abs_tol) - double(rel_tol) * std::abs(val);


    return tolerance_diff <= 0;
}


// abs is zero by default
template<typename V, typename R>
auto LeTolerant(V &&val, R &&rel_tol) {
    return LeTolerant(std::forward<V>(val), std::forward<R>(rel_tol), 0.0);
}


