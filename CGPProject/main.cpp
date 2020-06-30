#include <iostream>
#include "my_fourier.h"
#include "mycgp.h"


int main() {
    std::cout << " Started 1 " << std::endl;

    MyFourierClass f; 
    f.load_from_csv("newwaves.csv");
    
    f.forward_fft(10001);
    f.get_prev_output();
    f.inverse_fft(10001);
    std::cout << " Started 2 " << std::endl;
    //cgpWrapper::my_runCGP();
}
