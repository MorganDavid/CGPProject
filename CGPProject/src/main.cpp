#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"
#include "fftw3.h"

int main() {
    dataSet* ds = initialiseDataSetFromFile("datasets/1000withphase.csv");
    MyFourierClass f(100, ds);
    int terms = 3; 
    int* terms_arr = new int[terms] {1, 2, 5};
    f.set_terms_array(terms, terms_arr);
    f.execute_extract_harmonics(terms);
    
    f.write_harmonics_to_csv("harmonics");
    f.write_to_csv_1d("amps.csv",&f.get_amplitude_list()[0],f.get_amplitude_list().size());
    f.write_to_csv_1d("freq.csv",&f.get_frequency_list()[0],f.get_frequency_list().size());
    auto ya = f.getSynthesisWithHarmonics(terms);
    f.write_to_csv_1d("syntehsiz.csv", &ya[0],ya.size());

    // cgpWrapper::initializeParams();
    //cgpWrapper::harmonic_runCGP_wave("datasets/pjme_2_years.csv");
    //cgpWrapper::my_runCGP(initialiseDataSetFromFile("complex-300pnts.csv"));
}
