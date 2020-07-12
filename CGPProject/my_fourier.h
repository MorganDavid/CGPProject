#include "fftw3.h"
#include "fstream"
#include <complex>

template<class M>
struct myMatrix {// remeber to set height width and size. 
    M** data=nullptr;
    size_t height=0;
    size_t width=0;
};

class MyFourierClass {
    myMatrix<double> dataset; // Stores the matrix to be Fourier transformed. TODO: change this to use pointer matrix not boost.
    double Fs; //  sample rate 
public:
    MyFourierClass(double Fs) : Fs{ Fs } {};
    virtual ~MyFourierClass() {};

    void forward_fft(const int bins, const size_t L, fftw_complex* out) const;
    void fourier_series(const fftw_complex* freq_spect, const int terms, const size_t L, double** out_sin, double** out_cos) const;
    void synthesise_from_waves(const int terms, const size_t L,const double* const * sin_mat, const double*const* cos_mat, double** out_synthesis) const; //TODO: learn how to make matricies const properly.

    void inverse_fft(const int terms, fftw_complex* input, double* output) const;

    void load_from_csv(const std::string file_dir);

    static void extract_harmonics(std::string file_dir, double Fs, int terms, double** synthesized);

    //Template functions have to be inline (i think)
    template<class T>
    inline void write_to_csv(const char* file_dir, T** arr, const int width, const int height) const {
        std::ofstream out(file_dir);

        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++)
                out << arr[r][c] << ',';
            out << '\n';
        }
        out.close();
    };   
    /*
    Csv writes for 1d arrays
    */
    //Template functions have to be inline (i think)
    template<class T>
    inline void write_to_csv_1d(const char* file_dir, const T * const arr, const int L) const {
        std::ofstream out(file_dir);

        for (int r = 0; r < L; r++) {
            out << arr[r] << ',' << '\n';
        }
        out.close();
    };
    //TODO: make a new implementation of this to generalize to any type with tostring ability. 
    inline void write_to_csv(const char* file_dir, fftw_complex* arr, const int L) const {
        std::ofstream out(file_dir);

        for (int r = 0; r < L; r++) {
           out << arr[r][0] << ":" << arr[r][1] << ',' << "\n";
        }
        out.close();
    };

};
