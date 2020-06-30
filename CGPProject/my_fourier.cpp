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


/// <summary>
/// Computes forward fourier transform on this->dataset
/// </summary>
/// <returns>Frequency spectrum as fftw_complex* is put into this->prev_output</returns>
void MyFourierClass::forward_fft(int bins) {
    // TODO: test if dataset is intitialized first.
    // TEMPORARY: get first dimension from this for testing. 
    double* col = myhelpers::getColFromMatrix(dataset, 0);
    
    int num_rows = 10001;

    this->prev_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows);

    fftw_plan p;
    p = fftw_plan_dft_r2c_1d(num_rows, col, this->prev_output, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);

};

/// <summary>
/// Calculates the fourier series representation of a frequency spectrum using the largest "terms" amplitudes in the spectrum.
/// TODO: maybe look into using linear algebra library to speed this up. 
/// </summary>
/// <param name="freq_spect">output from fft</param>
/// <param name="terms">number of terms to reconstruct with.</param>
void MyFourierClass::fourier_series(fftw_complex* freq_spect, int terms) {
    using namespace std::complex_literals;
   /*
    double Fs = 10001; // sampel rate Hz
    double L = 1000; // # of samples.
    

    double* w;
    w = new double[L];
    for (int i = 0; i < L; ++i) {
        double new_time = i * (Fs / L); //<- Don't understand why I have to do this. 
        w[i] = 2 * M_PI * new_time;
    }
    */
    // Sort in order of amplitude
   /* [out, idx] = sort(yf, 'descend');
    yf_sorted = out;
    w = w(idx);*/
    /*
    double* out;
    out = new double[L];

    std::complex<double>* comp_freq_spect = reinterpret_cast<std::complex<double>*>(freq_spect); // convert to std::complex
    for (int k = 1; k < terms; ++k) {
        out[k] = out[k] + comp_freq_spect[k] * exp(1i * w[k] ); // *t here
    }
    */
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

// Load the dataset to be Fourier transfomred into this->dataset.
void MyFourierClass::load_from_csv(std::string file_dir) {
    const int num_cols = 1;
    const int num_rows = 10001;
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
            ss >> this->dataset(r,c);
        }

    }
};
#endif