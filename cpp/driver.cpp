#include "typedefs.h"

#include "cilk.h"
#include "HermiteRunner.h"
#include "Naive.h"
#include <argparse/argparse.hpp>
#include <experimental/mdspan>
#include <iostream>
#include <fstream>
#include <optional>

// TODO
using namespace ahr;

int main(int argc, const char *argv[]) {
    std::cout << "Hello!" << std::endl;

    argparse::ArgumentParser arguments("ahr");

    arguments.add_argument("X")
            .help("Size of X, Y for domain")
            .scan<'i', Dim>()
            .default_value(Dim{256});

    arguments.add_argument("M")
            .help("Total number of moments")
            .scan<'i', Dim>()
            .default_value(Dim{100});

    arguments.add_argument("nr")
            .help("Number of timesteps")
            .scan<'i', Dim>()
            .default_value(Dim{20});

    arguments.add_argument("inFile")
            .help("The path to the file that contains initial moment data")
            .default_value("");

    try {
        arguments.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << arguments;
        return 1;
    }

    auto X = arguments.get<Dim>("X"),
            M = arguments.get<Dim>("M"),
            nr = arguments.get<Dim>("nr");

    std::optional<std::ifstream> inFile{};
    if (auto fPath = arguments.get("inFile"); !fPath.empty()) {
        inFile = std::ifstream{fPath};
    }

    Naive naive{std::cout, M, X, X};
    HermiteRunner &runner = naive;

    stdex::mdarray<Real, stdex::dextents<size_t, 3u>> initialMoments{M, X, X};
    for (int m = 0; m < M; ++m) {
        for (int x = 0; x < X; ++x) {
            for (int y = 0; y < X; ++y) {
                if (inFile)
                    *inFile >> initialMoments(m, x, y);
                else
                    initialMoments(m, x, y) = 0;
            }
        }
    }

    runner.init(initialMoments, nr, 0.01);
    runner.run();
    auto moments = runner.getFinalValues();
}
