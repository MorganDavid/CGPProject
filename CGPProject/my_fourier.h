#include "fftw3.h"
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class MyFourierClass {
    matrix<double> dataset; // Stores the matrix to be Fourier transformed.
    fftw_complex* prev_output; // The answer to the last completed fourier calculation.
    double* inverse_output;
public:
    virtual ~MyFourierClass() {
        
        std::cout << " destroyed " << std::endl;
    };

    void forward_fft(int);
    void findTopKHarms(fftw_complex* freq_spect, int L, int terms, int* out_inds);
    void fourier_series(fftw_complex* freq_spect, int terms, int L, double** out_sin, double** out_cos);
    void inverse_fft(int);
    
    void load_from_csv(std::string file_dir);

    inline fftw_complex* get_prev_output() {
        return(this->prev_output);
    }
};
