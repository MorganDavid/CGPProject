#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"
#include "fftw3.h"

int main() {
   
    /*MyFourierClass f(300,"complex-300pnts.csv");

    f.execute_extract_harmonics(5);

    f.write_harmonics_to_csv("harmonics");
    f.write_to_csv_1d("amps",&f.get_amplitude_list()[0],f.get_amplitude_list().size());
    */
    cgpWrapper::initializeParams();
    cgpWrapper::harmonic_runCGP("complex-300pnts.csv");
   // cgpWrapper::my_runCGP(initialiseDataSetFromFile("ez_mode_sinwaves.csv"));
}
