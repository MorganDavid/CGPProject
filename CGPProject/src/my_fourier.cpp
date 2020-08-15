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
    // TODO: test if dataset is intitialized first.
    // TEMPORARY: get first dimension from this for testing.
    std::vector<double> col = myhelpers::getColFromMatrix(dataset, 0);
    MyFourierClass::write_to_csv_1d("blah.csv", &col[0],L);
    fftw_plan p;
    p = fftw_plan_dft_r2c_1d(L, &col[0], out, FFTW_ESTIMATE);

    fftw_execute(p);
    //divide output by L
    for (int i = 0; i < L / 2; i++) {
        out[i][0] = 2*abs(out[i][0]) / L;
        out[i][1] = 2*abs(out[i][1]) / L;
    }
    MyFourierClass::write_to_csv("out.csv", out, L);

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
void MyFourierClass::fourier_series(const std::vector<std::complex<double>> freq_spect_cmplx, const std::vector<double> amplitudeList, const int terms, const double Fs, const size_t L, double** out_sin, double** out_cos) {
    // -- Find strongest k amplitudes in freq_spect --
    std::vector<double> amp(amplitudeList);
    //write_to_csv_1d("amp_spec.csv", &amp[0], (int)L/2);

    std::vector<int> harms_idx = myhelpers::maxk(amp, terms); // indicies of top harmonics

    //std::vector<double> t(L * 2); // timestep vector TODO: remove this. (leaving for reference later.
    //std::generate(t.begin(), t.end(), [Fs]() { static int i = 0;  return ((1.0/Fs) * (double)i++ ); });

    // -- Extract component waves --
    for (int i = harms_idx.size() - 1; i >= 0; i--) {
        int k = harms_idx[i];
        double pr = 2 * M_PI * (k * (double)(Fs / L));//input, assuming regular sampling.

        out_cos[i] = new double[L];
        out_sin[i] = new double[L];

        for (int x = 0; x < L; x++) {
            double t = 1.0 / Fs * x;
            double cos_ans = (double)((freq_spect_cmplx[k]).real() * cos(pr * t));  //real(2*yf(k)/L)*cos(w*t(x));
            double sin_ans = (double)((freq_spect_cmplx[k]).imag() * sin(pr * t));
            out_cos[i][x] = cos_ans;
            out_sin[i][x] = sin_ans;
        }
    }
}

/// <summary>
/// Uses the results from `fourier_series` to synthesise the original function using numbers of harmonics UP TO `terms`. 
/// </summary>
/// <param name="out_synthesis">Must be initialized with first dimension size `terms`. Indexing out_synthesis[0] will give you the result of the fourier series with 1 harmonic,
/// out_synthesis[1] will give you the result with 2 harmonics and so on up to `terms`.</param>
void MyFourierClass::synthesise_from_waves(const int terms, const size_t L,const double*const* sin_mat, const double* const * cos_mat, double** out_synthesis) {
    for (int t = 1; t < terms+1; t++){
        out_synthesis[t-1] = new double[L]();
        for (int i = 0; i < t; i++) {
            for (int z = 0; z < L; z++) {
                double ans = out_synthesis[t-1][z] + sin_mat[i][z] + cos_mat[i][z];
                out_synthesis[t-1][z] = ans;
            }
        }
    }
    
}

std::vector<double> MyFourierClass::calculate_amplitude_list(std::vector<std::complex<double>> freq_spect) {
    std::vector<double> amp(freq_spect.size()); // amplitude vector
    for (int i = 0; i < freq_spect.size(); i++) {
        double abbed = std::abs(freq_spect[i]);
        amp[i] = (abbed); // = 2.*|freq_spec(1:L/2)./L|
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
    this->out_cos = new double* [terms];
    this->out_sin = new double* [terms];

    execute_forward_fft(L);
    this->frequencyList = myhelpers::fftw_complex2std_complex(this->freq_spect, L / 2);
    // Fix frequencies to be correct (dont understand this problem).
   
    this->amplitudeList = calculate_amplitude_list(frequencyList);
    //TODO: add phase here
    fourier_series(this->frequencyList, this->amplitudeList, terms, this->Fs, L, out_sin, out_cos);

    this->execute_synthesise_from_waves(terms, out_sin, out_cos);
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