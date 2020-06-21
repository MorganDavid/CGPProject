#include <iostream>
#include "fftw3.h"

void fourier_test() {
    {
        fftw_complex* in, * out;
        fftw_plan p;
        
        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 10);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 10);
        p = fftw_plan_dft_1d(10, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
       fftw_execute(p); /* repeat as needed */
        
        fftw_destroy_plan(p);
        fftw_free(in); fftw_free(out);
    }
}