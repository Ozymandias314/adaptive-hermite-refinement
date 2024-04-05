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

    arguments.add_argument("save-interval")
            .help("How often the code should write out the results")
            .scan<'i', Dim>()
            .default_value(Dim{1000});

    try {
        arguments.parse_args(argc, argv);
        if (arguments.get<Dim>("M") < 4) throw std::invalid_argument("At least 4 moments required");
    }
    catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << arguments;
        return 1;
    }

    auto X = arguments.get<Dim>("X"),
            M = arguments.get<Dim>("M"),
            nr = arguments.get<Dim>("nr"),
            saveInterval = arguments.get<Dim>("save-interval");

    Naive naive{std::cout, M, X, X};
    HermiteRunner &runner = naive;

    runner.init(nr);
    runner.run(saveInterval);
    auto aPar = runner.getFinalAPar();

    for (int x = 0; x < X; ++x) {
        for (int y = 0; y < X; ++y) {
            std::cout << aPar(x, y) << " ";
        }
        std::cout << std::endl;
    }
}
