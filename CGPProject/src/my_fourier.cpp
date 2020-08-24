#ifndef MY_FOURIER
#define MY_FOURIER
#define _USE_MATH_DEFINES

#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "my_fourier.h"
#include <cmath>
#include <complex>
#include "my_struct_definitions.c"
#include "myhelpers.h"

/// <summary>
/// Computes forward fourier transform on this->dataset
/// </summary>
/// <returns>Frequency spectrum as fftw_complex* is put into this->prev_output</returns>
void MyFourierClass::forward_fft(const int bins, const size_t L, myMatrix<double> dataset, fftw_complex* out) {
    // Only apply foui
    std::vector<double> col = myhelpers::getColFromMatrix(dataset, 0);
    MyFourierClass::write_to_csv_1d("blah.csv", &col[0], L);
    fftw_plan p;
    p = fftw_plan_dft_r2c_1d(L, &col[0], out, FFTW_ESTIMATE);

    fftw_execute(p);
    MyFourierClass::write_to_csv("before.csv", out, L/2);

    // divide output by L
    for (int i = 0; i < L / 2; i++) {
        out[i][0] = ((out[i][0]) / L);
        out[i][1] = ((out[i][1]) / L);
    }
    MyFourierClass::write_to_csv("out.csv", out, L/2);

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
void MyFourierClass::fourier_series(const std::vector<std::complex<double>> freq_spect_cmplx, const std::vector<double> amplitudeList, const int terms, const double Fs, const size_t L, std::complex<double>* out_synthesis) {
    using namespace std::complex_literals;
    // -- Find strongest k amplitudes in freq_spect --
    std::vector<double> amp(amplitudeList);

    std::vector<int> harms_idx = myhelpers::maxk(amp, terms); // indicies of top harmonics

    for (int x = 0; x < harms_idx.size(); x++) {
        int k = harms_idx[x];
        double pr = 2 * M_PI * (k * (double)(Fs / L)); // input, assuming regular sampling.
        
        for (int i = 0; i < L; i++) {
            double t = 1.0 / Fs * i;
            out_synthesis[i] = out_synthesis[i] + freq_spect_cmplx[k] * exp(1i * pr * t);
        }
    }
}

std::vector<double> MyFourierClass::calculate_amplitude_list(std::vector<std::complex<double>> freq_spect) {
    std::vector<double> amp(freq_spect.size()); // amplitude vector
    for (int i = 0; i < freq_spect.size(); i++) {
        double abbed = std::abs(freq_spect[i]);
        amp[i] = (abbed); 
    }
    return amp; 
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
// DISCLAIMER: copied some of this function's code from cgp.c
// First line of csv must be inputs,outputs,samples. Same as cgp library.
// Also allocates data for this->dataset.
// TODO: this function has stopped working, fix it. 
void MyFourierClass::load_from_csv(const std::string file_dir) {
    std::ifstream file(file_dir);
    std::string line;
    std::getline(file, line);
    int num_inputs = -1;
    int num_outputs = -1;
    int num_rows= -1;//not using size_t bc it doesn't work with sscanf_s
    sscanf_s(line.c_str(), "%d,%d,%d", &num_inputs, &num_outputs, &num_rows);
    int num_cols = num_inputs + num_outputs;
    this->dataset.data = new double* [num_rows];
    this->dataset.height = num_rows;

    for (int r = 0; r < num_rows; r++)
    {
        std::string line;
        std::getline(file, line);

        if (!file.good()) {
            std::cout<< ("ERROR "+file_dir + " at row " + std::to_string(r))<<std::endl;
            break;
        }

        std::stringstream iss(line);
        this->dataset.data[r] = new double [num_cols];
        for (int c = 0; c < num_cols; c++)
        {
            std::string v;
            std::getline(iss, v, ',');
            if (!iss.good()) {
                std::cout << ("ERROR" + file_dir + " at " + std::to_string(r) + "," + std::to_string(c))<<std::endl;
                break;
            }
            std::stringstream ss(v);
            ss >> this->dataset.data[r][c];
        }

    }
    this->dataset.width = num_cols;
}

// Converts from the cgp library dataset to myMatrix using copy. 
void MyFourierClass::setDatasetFromCGP(struct dataSet* data)
{
    this->dataset.height = data->numSamples;
    this->dataset.width = data->numOutputs;
    this->dataset.data = new double*[this->dataset.height];
    for (int i = 0; i < this->dataset.height; i++) {
        this->dataset.data[i] = new double [this->dataset.width];
        for (int x = 0; x < this->dataset.width; x++) {

            this->dataset.data[i][x]=data->outputData[i][x];
        }
    }
};

/// <summary>
/// Allows the getharmonic function to be run to get harmonics from the dataset in this object. 
/// </summary>
/// <param name="file_dir"></param>
/// <param name="Fs"></param>
/// <param name="terms"></param>
/// <param name="out_synth"> allocate with  `double** out_synth = new double* [terms];`</param>
void MyFourierClass::execute_extract_harmonics(int terms)
{
    int L = this->dataset.height;
    this->num_of_terms = terms;

    execute_forward_fft(L);
    this->frequencyList = myhelpers::fftw_complex2std_complex(this->freq_spect, L / 2);
    // Fix frequencies to be correct (dont understand this problem).
   
    this->amplitudeList = calculate_amplitude_list(frequencyList);

    this->harmonic_output.data = new double* [terms];
    this->harmonic_output.height = terms;
    this->harmonic_output.width = this->dataset.height;

    // Each row of `harmonic_output`.data will represent an incrementing number of harmonics.
    for (int i = 0; i < terms; i++) {
        int run_terms = i;
        if (terms_array != NULL) run_terms = this->get_terms_array()[i]; // if terms doesn't always increment by one. 

        auto out_synthesis = new std::complex<double> [L];
        fourier_series(this->frequencyList, this->amplitudeList, run_terms, this->Fs, L, out_synthesis);
        
        this->harmonic_output.data[i] = new double [L];
        for (int x = 0; x < L; x++) {
            this->harmonic_output.data[i][x] = out_synthesis[x].real();
        }
    }
    write_to_csv("harmonics_output_data.csv", this->harmonic_output.data, harmonic_output.width, harmonic_output.height);
}

/// <summary>
/// Only run this after `execute_extract_harmonics`. Gets synthesis with `num_harmonics` from the harmonics_output.
/// </summary>
/// <param name="num_harmonics"></param>
/// <returns></returns>
std::vector<double> MyFourierClass::getSynthesisWithHarmonics(int num_harmonics)
{
    std::vector<double> out(harmonic_output.data[num_harmonics - 1], harmonic_output.data[num_harmonics - 1]+harmonic_output.width);
    return out;
};

#endif