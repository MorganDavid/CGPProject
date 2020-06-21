#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "myhelpers.cpp";

class MyFourierClass {
public:
    virtual ~MyFourierClass() {
        fftw_free(prev_output);
        delete[](dataset);
        std::cout << " destroyed " << std::endl;
    }

    //Computes forward fourier transform on this->dataset
    fftw_complex* forward_fft() {
        if (dataset == 0) { throw "Dataset in MyFourierClass must be intitilazed first."; }
        //TEMPORARY: get first dimension from this for testing. 
        double* col = myhelpers::getColFromMatrix(dataset,0,10001);

        int num_rows = 10001;

        this->prev_output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * num_rows);

        fftw_plan p;
        p = fftw_plan_dft_r2c_1d(num_rows, col, prev_output, FFTW_ESTIMATE);

        fftw_execute(p);

        fftw_destroy_plan(p);

    }

    // Load the dataset to be Fourier transfomred into this->dataset.
    void load_from_csv(std::string file_dir) {
        const int num_cols = 1;
        const int num_rows = 10001;
        this->dataset = new double* [num_rows]; //TODO: get size from file (put a header at the top of the file).

        std::ifstream file(file_dir);

        for (int r = 0; r < num_rows; r++)
        {
            std::string line;
            std::getline(file, line);

            if (!file.good()) {
                std::cout << "bad file (row): " << std::endl;
                break;
            }

            dataset[r] = new double[num_cols];

            std::stringstream iss(line);

            for (int c = 0; c < num_cols; c++)
            {
                std::string v;
                std::getline(iss, v, ',');
                if (!iss.good()) {
                    std::cout << "bad file (col): " << file_dir << std::endl;
                    break;
                }
                std::stringstream ss(v);
                ss >> dataset[r][c];
            }

        }
    }
    //This performs on complex numbers. 
   /* void fourier_test() {
        fftw_complex* in, *;
        fftw_plan p;

        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 10);

        this->prev_output= (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 10);

        p = fftw_plan_dft_1d(10, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        fftw_execute(p); 

        //print out array
        for (int i = 0; i < sizeof(out); i++) {
            std::cout << out[i] << ": " << std::endl;
        }

        fftw_destroy_plan(p);
        fftw_free(in);

    }
    */
    /*
    The results of the previous Fourier calculation. 
    */
    fftw_complex* get_prev_output() {
        return(this->prev_output);
    }

/*
   PRIVATE Implementations
*/
private:
    double** dataset = 0; // Stores the matrix to be Fourier transformed.
    fftw_complex* prev_output = 0; // The answer to the last completed fourier calculation.

  


};