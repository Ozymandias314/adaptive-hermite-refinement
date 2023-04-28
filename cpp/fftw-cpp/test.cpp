#include <complex.h>
#include <fftw3.h>
#include <numbers>
#include <iostream>

void print(const fftw_complex *arr, int N) {
    for (int i = 0; i < N; ++i) {
        std::cout << "("<<arr[i][0] << "," << arr[i][1] << ") ";
    }
    std::cout << std::endl;
};

int main() {
    fftw_complex *in, *out, *out2;
    fftw_plan p;
    int N = 8;

     in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
    out2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pInv = fftw_plan_dft_1d(N, out, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int j = 0; j < N; ++j) {
        in[j][0] = std::cos(2.0 * std::numbers::pi * double(j) / double(N));
        in[j][1] = std::cos(2.0 * std::numbers::pi * (j+N/2) / N );
    }

    fftw_execute(p); /* repeat as needed */
    fftw_execute(pInv); /* repeat as needed */

    print(in, N);
    print(out, N);
    print(out2, N);

    auto *out2_cpp = reinterpret_cast<std::complex<double> *>(out2);
    for (int i = 0; i < N; ++i) {
        std::cout << out2_cpp[i] / double(N) << " ";
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(pInv);

    fftw_free(in);
    fftw_free(out);
    fftw_free(out2);
}