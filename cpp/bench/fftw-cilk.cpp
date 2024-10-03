#include "fftw-cpp/fftw-cpp.h"

#include <benchmark/benchmark.h>
#include <cilk/cilk.h>
#include <cstring>
#include <iostream>

static void copy(fftw::buffer const &in, fftw::buffer &out) {
  std::memcpy(out.data(), in.data(),
              in.size() * sizeof(fftw::buffer::element_type));
}

static void BM_copy_1D(benchmark::State &state) {
  size_t const X = state.range(0);

  fftw::buffer in{X};
  fftw::buffer out{X};

  // warmup
  copy(in, out);

  for (auto _ : state) {
    copy(in, out);
  }
}

static void BM_FFTW_1D(benchmark::State &state) {
  size_t const X = state.range(0);

  fftw::buffer in{X};
  fftw::buffer out{X};

  fftw_plan_with_nthreads(int(X));
  auto const plan = fftw::plan<>::dft(in, out, fftw::FORWARD, fftw::MEASURE);

  // warmup
  plan();

  for (auto _ : state) {
    plan();
  }
}

static void BM_copy_2D(benchmark::State &state) {
  size_t const X = state.range(0);
  size_t const Y = state.range(1);

  fftw::mdbuffer<2> in{X, Y};
  fftw::mdbuffer<2> out{X, Y};

  // warmup
  copy(in.container(), out.container());

  for (auto _ : state) {
    copy(in.container(), out.container());
  }
}

static void BM_FFTW_2D(benchmark::State &state) {
  size_t const X = state.range(0);
  size_t const Y = state.range(1);

  fftw::mdbuffer<2> in{X, Y};
  fftw::mdbuffer<2> out{X, Y};


  fftw_plan_with_nthreads(int(X*Y));
  auto const plan = fftw::plan<2>::dft(in, out, fftw::FORWARD, fftw::MEASURE);

  // warmup
  plan();

  for (auto _ : state) {
    plan();
  }
}

auto constexpr K = 1024;
auto constexpr M = K * K;

BENCHMARK(BM_copy_1D)->Range(16 * M, 128 * M)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_FFTW_1D)->Range(16 * M, 128 * M)->Unit(benchmark::kMillisecond);

BENCHMARK(BM_copy_2D)
    ->ArgsProduct({{4 * K, 8 * K, 16 * K}, {4 * K, 8 * K, 16 * K}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_FFTW_2D)
    ->ArgsProduct({{4 * K, 8 * K, 16 * K}, {4 * K, 8 * K, 16 * K}})
    ->Unit(benchmark::kMillisecond);

// TODO pthreads, openmp, cilk -> compare all
#ifndef PTHREADS
void parallel_for(void *(*work)(char *), char *jobdata, size_t elsize,
                  int njobs, void *data) {
  // std::cout << "parallel_for: " << njobs << std::endl;
  cilk_for(int i = 0; i < njobs; ++i) { work(jobdata + i * elsize); }
}
#endif

int main(int argc, char **argv) {
  fftw_init_threads();
  fftw_plan_with_nthreads(K);
  fftw_threads_set_callback(parallel_for, nullptr);

  // make sure cilk isn't initialized and deinitialized for each benchmark
  cilk_scope {
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
  }
}
