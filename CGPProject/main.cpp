#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"
#include "fftw3.h"

int main() {
  /*  MyFourierClass f(300,"complex-300pnts.csv");

    f.execute_extract_harmonics(5);

    f.write_harmonics_to_csv("harmonics");

    const int terms = 8;
    const int L = 300;
    double** out_cos = new double* [terms];
    double** out_sin = new double* [terms];
    double** out_synth = new double* [terms];
    fftw_complex* freq_spec = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * L/2);

    MyFourierClass f(300); 
    f.load_from_csv("complex-300pnts.csv");

    f.forward_fft(L, L, freq_spec);

    f.fourier_series(freq_spec, terms, L, out_sin, out_cos);
    f.write_to_csv("out_sin.csv", out_sin, L, terms);
    f.write_to_csv("out_cos.csv", out_cos, L, terms);

   // f.inverse_fft(L, freq_spec);
    f.synthesise_from_waves(terms, L, out_sin, out_cos, out_synth);

    // Write to csv for debugging. 

    f.write_to_csv_1d("freq_spec.csv", freq_spec, L/2);

    f.write_to_csv("out_synth.csv", out_synth, L, terms);

    //Cleanup
    for (int i = 0; i < terms; i++) {
        delete out_cos[i];
        delete out_sin[i];
    //    delete out_synth[i];
    }

    delete[] out_cos;
    delete[] out_sin;*/
    //delete[] out_synth;



    cgpWrapper::my_runCGP();
}
