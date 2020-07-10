#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"
#include "fftw3.h"

int main() {
    const int terms = 8;
    const int L = 300;
    double** out_cos = new double* [terms];
    double** out_sin = new double* [terms];
    double** out_synth = new double* [terms];
    fftw_complex* freq_spec = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * L/2);

    MyFourierClass f; 
    f.load_from_csv("complex-300pnts.csv"); //TODO: this csv is in the cgp library format. do something to mak eit work with my code. 

    f.forward_fft(L, L, freq_spec);

    f.fourier_series(freq_spec, terms, L, out_sin, out_cos);
    //f.inverse_fft(L, freq_spec);
   // f.synthesise_from_waves(5, L, out_sin, out_cos, out_synth);

    // Write to csv for debugging. 
    f.write_to_csv<double>("out_sin.csv", out_sin, L, terms);
    f.write_to_csv<double>("out_cos.csv", out_cos, L, terms);

    f.write_to_csv("freq_spec.csv", freq_spec, L);
    
    //Cleanup
    for (int i = 0; i < terms; i++) {
        delete out_cos[i];
        delete out_sin[i];
      //  delete out_synth[i];
    }

    delete[] out_cos;
    delete[] out_sin;
    //delete[] out_synth;


    //cgpWrapper::my_runCGP();
}
