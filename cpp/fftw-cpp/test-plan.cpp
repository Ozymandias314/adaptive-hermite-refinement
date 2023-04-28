#include <fftw-cpp/fftw-cpp.h>

int main(){
    fftw::buffer buf1{100};
    auto plan1 = fftw::plan<1u>::dft(buf1, buf1, fftw::FORWARD, fftw::ESTIMATE);
    auto plan2 = std::move(plan1);

    plan1 = std::move(plan2);

    auto buf2 = std::move(buf1);
    buf1 = std::move(buf2);
}