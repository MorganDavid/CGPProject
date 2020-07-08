#include "fftw3.h"
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

class MyFourierClass {
    matrix<double> dataset; // Stores the matrix to be Fourier transformed.
    fftw_complex* prev_output; // The answer to the last completed fourier calculation.
    double* inverse_output;
    double Fs = 3000;
public:
    virtual ~MyFourierClass() {
        
        std::cout << " destroyed " << std::endl;
    };

    void forward_fft(int);
    void inverse_fft(int);
    void fourier_series(const fftw_complex* freq_spect, const int terms, const size_t L, double** out_sin, double** out_cos) const;
    
    void load_from_csv(std::string file_dir);

    inline fftw_complex* get_prev_output() {
        return(this->prev_output);
    }
    template<class T>
    void write_to_csv(const char* file_dir, T** file, const int width, const int height) const;
    template<>
    void write_to_csv(const char* file_dir, fftw_complex** file, const int width, const int height) const;

};
