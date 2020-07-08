#ifndef MY_FOURIER
#define MY_FOURIER
#define _USE_MATH_DEFINES

#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "my_fourier.h"
#include "myhelpers.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <complex>

/// <summary>
/// Computes forward fourier transform on this->dataset
/// </summary>
/// <returns>Frequency spectrum as fftw_complex* is put into this->prev_output</returns>
void MyFourierClass::forward_fft(int bins) {
    // TODO: test if dataset is intitialized first.
    // TEMPORARY: get first dimension from this for testing. 
    std::vector<double> col = myhelpers::getColFromMatrix(dataset, 0);

    int num_rows = 6000;

    this->prev_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows);

    fftw_plan p;
    p = fftw_plan_dft_r2c_1d(num_rows, &col[0], this->prev_output, FFTW_ESTIMATE);

    fftw_execute(p);

    //for (int i = 0; i < 10; i++) std::cout << prev_output[i][0] << "+"<< prev_output[i][1]<< std::endl;
    fftw_destroy_plan(p);

};

/// <summary>
/// Calculates the fourier series representation of a frequency spectrum using the largest "terms" amplitudes in the spectrum.
/// Define outputs like this: 
/// <code>
/// double** out_cos = new double* [terms];
/// double** out_sin = new double* [terms];
/// Cleanup not handled by this class !
/// </code>
/// </summary>
/// <param name="freq_spect">output from fft</param>
/// <param name="terms">number of terms to reconstruct with.</param>
/// TODO: I think this code produces the output flipped/backwards for some reason. 
void MyFourierClass::fourier_series(const fftw_complex* freq_spect, const int terms, const size_t L, double** out_sin, double** out_cos) const {
    // -- Find strongest k amplitudes in freq_spect --
    std::vector<std::complex<double>> freq_spect_cmplx = myhelpers::fftw_complex2std_complex(freq_spect, L);

    std::vector<double> amp(L); // amplitude vector
    for (int i = 0; i < L; i++) {
        double abbed = std::abs(freq_spect_cmplx[i]);
        amp[i] = 2 * abbed / (L * 2.0); // = 2.*|freq_spec(1:L/2)./L|
    }

    std::vector<int> harms_idx = myhelpers::maxk(amp, terms); // indicies of top harmonics

    
    //std::vector<double> t(L * 2); // timestep vector TODO: remove this. (leaving for reference later.
    //std::generate(t.begin(), t.end(), [Fs]() { static int i = 0;  return ((1.0/Fs) * (double)i++ ); });
    
    // -- Extract component waves --
    for (int i = 0; i < harms_idx.size(); i++) {
        int k = harms_idx[i];
        double pr = 2 * M_PI * ((k - 1.0) * (this->Fs / L*2.0));

        out_cos[i] = new double[L*2];
        out_sin[i] = new double[L*2];
        for (int x = 0; x < L*2; x++) {
            double t = 1.0 / this->Fs * x;
            out_cos[i][x] = (double)((2 * freq_spect_cmplx[k] / L*2).real() * cos(pr * t));  //real(2*yf(k)/L)*cos(w*t(x));
            out_sin[i][x] = (double)((2 * freq_spect_cmplx[k] / L*2).imag() * sin(pr * t));
        }
    }

}

/// <summary>
/// Uses the value in inverse_output to compute the inverse fft. result is stored in dataset. 
/// notes: data in inverse_output will be destroyed (dangling pointer)
/// </summary>
/// <param name="terms">Number of terms to use in the approximation.</param>
void MyFourierClass::inverse_fft(int terms) {
    if (prev_output == 0) { throw "prev_output in MyFourierClass must be intitilazed first."; }

    this->inverse_output = (double*)fftw_malloc(sizeof(double) * terms);

    fftw_plan p;
    p = fftw_plan_dft_c2r_1d(terms, this->prev_output, this->inverse_output, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);
};

template<>//TODO: make a new implementation of this to generalize to any type with tostring ability. 
void MyFourierClass::write_to_csv(const char* file_dir, fftw_complex** arr, const int width, const int height) const{
    std::ofstream out(file_dir);

    for (int r = 0; r < height; r++){ 
        for (int c = 0; c < width; c++)
            out << arr[r][c] << ',';
        out << '\n';
    }
    out.close();
};

template<class T>
void MyFourierClass::write_to_csv(const char* file_dir, T** arr, const int width, const int height) const {
    std::ofstream out(file_dir);

    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++)
            out << arr[r][c] << ',';
        out << '\n';
    }
    out.close();
};

// Load the dataset to be Fourier transfomred into this->dataset.
void MyFourierClass::load_from_csv(std::string file_dir) {
    const int num_cols = 1;
    const int num_rows = 6000;
    this->dataset.resize(num_rows,num_cols); //TODO: get size from file (put a header at the top of thecsv  file).

    std::ifstream file(file_dir);

    for (int r = 0; r < num_rows; r++)
    {
        std::string line;
        std::getline(file, line);

        if (!file.good()) {
            myhelpers::handleError(0, file_dir + " at row " + std::to_string(r));
            break;
        }

        std::stringstream iss(line);

        for (int c = 0; c < num_cols; c++)
        {
            std::string v;
            std::getline(iss, v, ',');
            if (!iss.good()) {
                myhelpers::handleError(0, (file_dir + " at " + std::to_string(r) + "," + std::to_string(c)));
                break;
            }
            std::stringstream ss(v);
            ss >> this->dataset(r, c);
        }

    };

};
#endif