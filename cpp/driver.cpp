#include "typedefs.h"

#include "cilk.h"
#include "HermiteRunner.h"
#include "Naive.h"
#include <argparse/argparse.hpp>
#include <experimental/mdspan>
#include <iostream>

// TODO
using namespace ahr;

int main(int argc, const char *argv[]) {
    std::cout << "Hello!" << std::endl;

    argparse::ArgumentParser arguments("ahr");

    arguments.add_argument("K")
            .help("Size of Ky, Kx for domain")
            .scan<'i', Dim>()
            .default_value(Dim{256});

    arguments.add_argument("M")
            .help("Number of moments")
            .scan<'i', Dim>()
            .default_value(Dim{100});

    arguments.add_argument("nr")
            .help("Number of timesteps")
            .scan<'i', Dim>()
            .default_value(Dim{20});


    try {
        arguments.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << arguments;
        return 1;
    }

    auto K = arguments.get<Dim>("K"),
            M = arguments.get<Dim>("M"),
            nr = arguments.get<Dim>("nr");


    Naive naive{std::cout};
    HermiteRunner &runner = naive;

    // TODO
    naive.init({}, nr, 0.01);
    naive.run();
    auto moments = naive.getFinalValues();
}