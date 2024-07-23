#include "Naive.h"
#include "util.hpp"

#include <cnpy.h>

#include <filesystem>
#include <gtest/gtest.h>

class NaiveMoments : public ::testing::Test {
  protected:
    // Ignore output
    std::ostringstream out;

    void TearDown() override {
        if (HasFailure()) { std::cout << "Full output:" << out.str() << std::endl; }
    }
};

TEST_F(NaiveMoments, BasicOT01) {
    ahr::Naive naive{out, 10, 16, 16};

    naive.init("OT01");
    naive.run(10, 0); // no saving

    // Tests are run from inside the `test` directory
    auto aParRefNpy = cnpy::npy_load("./_test_data/a_par_16_16_10.npy");
    ASSERT_EQ(aParRefNpy.word_size, sizeof(double));
    std::span<size_t, 2> extents{aParRefNpy.shape.data(), 2};
    ahr::Naive::ViewXY view{aParRefNpy.data<double>(), extents};

    auto aPar = naive.getFinalAPar();
    ASSERT_THAT(aPar.to_mdspan(), MdspanElementsAllClose(view, 1e-10, 1e-10));
}
