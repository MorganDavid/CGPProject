#include "fftw3.h"
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class MyFourierClass {
    matrix<double> dataset; // Stores the matrix to be Fourier transformed.
    fftw_complex* prev_output = 0; // The answer to the last completed fourier calculation.
    double* inverse_output = 0;
public:
    virtual ~MyFourierClass() {
        fftw_free(prev_output);
        std::cout << " destroyed " << std::endl;
    };

    void forward_fft(int);
    void inverse_fft(int);

    void fourier_series(fftw_complex*, int);

    void load_from_csv(std::string file_dir);

    inline fftw_complex* get_prev_output() {
        return(this->prev_output);
    }
};
