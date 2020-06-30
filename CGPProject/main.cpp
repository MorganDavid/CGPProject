#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"

int main() {
    std::cout << " Started 1 " << std::endl;

    MyFourierClass f; 
    f.load_from_csv("6kfunc.csv");
    
    f.forward_fft(6000);
    int terms = 5;
    double** out_cos = new double* [terms];
    double** out_sin = new double* [terms];
    f.fourier_series(f.get_prev_output(), terms, 3000, out_sin, out_cos);
        
    //f.get_prev_output();
    //f.inverse_fft(10001);
    //cgpWrapper::my_runCGP();
}
