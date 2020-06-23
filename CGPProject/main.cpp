#include <iostream>
#include "my_fourier.cpp";
#include "mycgp.cpp"

int main() {
    std::cout << " Started 2 " << std::endl;
    MyFourierClass f; 
    f.load_from_csv("C:\\Users\\morga\\source\\repos\\CGPProject\\CGPProject\\newwaves.csv");
    f.forward_fft();
    f.get_prev_output();

    //cgpWrapper::my_runCGP();
}
