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

/*
 NOTE: All data which goes into this class must have EVEN LENGTH in the first dimension.
 L shoudl always be hte legnth of the original signal.
*/

/// <summary>
/// Computes forward fourier transform on this->dataset
/// </summary>
/// <returns>Frequency spectrum as fftw_complex* is put into this->prev_output</returns>
void MyFourierClass::forward_fft(const int bins, const size_t L, fftw_complex* out) const {
    // TODO: test if dataset is intitialized first.
    // TEMPORARY: get first dimension from this for testing. 
    std::vector<double> col = myhelpers::getColFromMatrix(dataset, 1);

    fftw_plan p;
    p = fftw_plan_dft_r2c_1d(L, &col[0], out, FFTW_ESTIMATE);

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
    std::vector<std::complex<double>> freq_spect_cmplx = myhelpers::fftw_complex2std_complex(freq_spect, L/2);

    std::vector<double> amp(L/2); // amplitude vector
    for (int i = 0; i < L/2; i++) {
        double abbed = std::abs(freq_spect_cmplx[i] / (L));
        amp[i] = 2 * abbed; // = 2.*|freq_spec(1:L/2)./L|
    }
    write_to_csv("amp_spec.csv", &amp[0], (int)L/2);

    std::vector<int> harms_idx = myhelpers::maxk(amp, terms); // indicies of top harmonics

    //std::vector<double> t(L * 2); // timestep vector TODO: remove this. (leaving for reference later.
    //std::generate(t.begin(), t.end(), [Fs]() { static int i = 0;  return ((1.0/Fs) * (double)i++ ); });

    // -- Extract component waves --
    for (int i = 0; i < harms_idx.size(); i++) {
        int k = harms_idx[i];
        double pr = 2 * M_PI * ((k - 1.0) * (this->Fs / L));

        out_cos[i] = new double[L];
        out_sin[i] = new double[L];

        for (int x = 0; x < L; x++) {
            double t = 1.0 / this->Fs * x;
            double cos_ans = (double)((2 * freq_spect_cmplx[k] / L).real() * cos(pr * t));  //real(2*yf(k)/L)*cos(w*t(x));
            double sin_ans = (double)((2 * freq_spect_cmplx[k] / L).imag() * sin(pr * t));
            out_cos[i][x] = cos_ans;
            out_sin[i][x] = sin_ans;
        }
    }
}
/// <summary>
/// 
/// </summary>
/// <param name="terms"></param>
/// <param name="L"></param>
/// <param name="sin_mat"></param>
/// <param name="cos_mat"></param>
/// <param name="out_synthesis">Must be initialized with first dimension size "terms". other dimension filled by function..</param>
void MyFourierClass::synthesise_from_waves(const int terms, const size_t L,const double*const* sin_mat, const double* const * cos_mat, double** out_synthesis) const {
    for (int t = 1; t < terms; t++){
        out_synthesis[t] = new double[L];
        for (int i = 0; i < t; i++) { 
            for (int z = 0; z < L; z++) {
                out_synthesis[t][z] = out_synthesis[t][z] + sin_mat[i][z] + cos_mat[i][z];
            }
        }
    }
}

/// <summary>
/// Uses the value in inverse_output to compute the inverse fft. result is stored in dataset. 
/// notes: data in input will be destroyed (leaves dangling pointer)
/// </summary>
/// <param name="terms">Number of terms to use in the approximation.</param>
void MyFourierClass::inverse_fft(const int terms, fftw_complex* input, double* output) const {
    output = (double*)fftw_malloc(sizeof(double) * terms);

    fftw_plan p;
    p = fftw_plan_dft_c2r_1d(terms, input, output, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);
};

// Load the dataset to be Fourier transfomred into this->dataset.
// DISCLAIMER: copied most of this code from cgp.c
// First line of csv must be inputs,outputs,samples. Same as cgp library.
void MyFourierClass::load_from_csv(const std::string file_dir) {
    std::ifstream file(file_dir);
    std::string line;
    std::getline(file, line);
    int num_inputs = -1;
    int num_outputs = -1;
    int num_rows= -1;//not using size_t bc it doesn't work with sscanf_s
    sscanf_s(line.c_str(), "%d,%d,%d", &num_inputs, &num_outputs, &num_rows);
    int num_cols = num_inputs + num_outputs;
    this->dataset.resize(num_rows, num_cols); 

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