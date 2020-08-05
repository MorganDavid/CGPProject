#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"
#include "fftw3.h"

int main() {
   
   /*  MyFourierClass f(500,"new_test.csv");
     
    f.execute_extract_harmonics(5);

    f.write_harmonics_to_csv("harmonics.csv");
    f.write_to_csv_1d("amps.csv",&f.get_amplitude_list()[0],f.get_amplitude_list().size());
    f.write_to_csv_1d("freq.csv",&f.get_frequency_list()[0],f.get_frequency_list().size());
  */
    cgpWrapper::initializeParams();
    cgpWrapper::harmonic_runCGP("datasets/1000pntsAt100Fs.csv");
    //cgpWrapper::my_runCGP(initialiseDataSetFromFile("complex-300pnts.csv"));
}
