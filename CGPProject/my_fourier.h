#include "fftw3.h"
#include "fstream"
#include <complex>
#include "dataSetDeclaration.h"

/*
 NOTE: All data which goes into this class must have EVEN LENGTH in the first dimension.
 L shoudl always be hte legnth of the original signal.
 Class object must be destroyed and reinstantiated every time you want to run on new data otherweise it will leak memory.
*/

template<class M>
class myMatrix {// remeber to set height width and size. 
public: 
    M** data=nullptr;
    size_t height=0;
    size_t width=0;
};

class MyFourierClass {
public:
    void execute_extract_harmonics(int terms);
    void getWithHarmonics(int number);

    void setDatasetFromCGP(struct dataSet* data);
    void load_from_csv(const std::string file_dir);

    MyFourierClass(double _Fs, dataSet* _dataset) :Fs{ _Fs } {
        this->setDatasetFromCGP(_dataset);
    };

    MyFourierClass(double _Fs, std::string file_dir) :Fs{ _Fs } {
        this->load_from_csv(file_dir);
    };

    virtual ~MyFourierClass() {
        if (dataset.height > 0) {
            for (int i = 0; i < dataset.height; i++) {
                delete dataset.data[i];
            }
            delete[] dataset.data;
        }
        // will memory leak if freq_spect has been set but not harmonic_output
        if (harmonic_output.height > 0) {
            for (int i = 0; i < harmonic_output.height; i++) {
                delete harmonic_output.data[i];
            }
            delete[] harmonic_output.data;
            fftw_free(freq_spect);
        }
    };

    inline void write_harmonics_to_csv(std::string file_dir) const {
        MyFourierClass::write_to_csv<double>(file_dir,this->harmonic_output.data, (int)this->harmonic_output.width, (int)this->harmonic_output.height);
    }

private:
    myMatrix<double> dataset; // Stores the matrix to be Fourier transformed.
    double Fs; //  sample rate 
    myMatrix<double> harmonic_output;
    fftw_complex* freq_spect;

    static void forward_fft(const int bins, const size_t L, myMatrix<double> dataset, fftw_complex* out);
    inline void execute_forward_fft(int bins) {
        size_t L = dataset.height;
        this->freq_spect = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * L / 2);
        forward_fft(bins, L, this->dataset, this->freq_spect);
    };

    static void fourier_series(const fftw_complex* freq_spect, const int terms, const double Fs, const size_t L, double** out_sin, double** out_cos);

    static void synthesise_from_waves(const int terms, const size_t L,const double* const * sin_mat, const double*const* cos_mat, double** out_synthesis);
    inline void execute_synthesise_from_waves(const int terms, const double* const* out_sin, const double* const* out_cos) {
        this->harmonic_output.data = new double* [terms];
        this->harmonic_output.height = terms;
        this->harmonic_output.width = this->dataset.height;
        synthesise_from_waves(terms, this->dataset.height, out_sin, out_cos, this->harmonic_output.data);
    };

    void inverse_fft(const int terms, fftw_complex* input, double* output) const;


    //Template functions have to be inline (i think)
    template<class T>
    static inline void write_to_csv(std::string file_dir, T** arr, const int width, const int height) {
        std::ofstream out(file_dir+".out");

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
    static inline void write_to_csv_1d(const char* file_dir, const T * const arr, const int L) {
        std::ofstream out(file_dir);

        for (int r = 0; r < L; r++) {
            out << arr[r] << ',' << '\n';
        }
        out.close();
    };

    //TODO: make a new implementation of this to generalize to any type with tostring ability. 
    static inline void write_to_csv(const char* file_dir, fftw_complex* arr, const int L) {
        std::ofstream out(file_dir);

        for (int r = 0; r < L; r++) {
           out << arr[r][0] << ":" << arr[r][1] << ',' << "\n";
        }
        out.close();
    };

};
