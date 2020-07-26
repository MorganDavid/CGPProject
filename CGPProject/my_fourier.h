#include "fftw3.h"
#include "fstream"
#include <complex>
#include <vector>
/*
 NOTE: All data which goes into this class must have EVEN LENGTH in the first dimension.
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
    std::vector<double> getSynthesisWithHarmonics(int num_harmonics);

    void setDatasetFromCGP(struct dataSet* data);
    void load_from_csv(const std::string file_dir);
    
    // These functions get the coefficients which represent the three parts of each wave.
    // up to the num_harmonics. your responsibility to pick sensible number for num_harmonics.
    std::vector<double> getAmplitudeList(int num_harmonics);
    std::vector<double> getFrequencyList(int num_harmonics);
    std::vector<double> getPhaseList(int num_harmonics);
    
    MyFourierClass(double _Fs, dataSet* _dataset) :Fs{ _Fs } { // TODO: this doesn't load the datset
        this->setDatasetFromCGP(_dataset);
        this->freq_spect = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dataset.height);

    };

    MyFourierClass(double _Fs, std::string file_dir) :Fs{ _Fs } {
        this->load_from_csv(file_dir);
        this->freq_spect = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dataset.height);
    };

    virtual ~MyFourierClass() {
       if (dataset.height > 0) {
            for (int i = 0; i < dataset.height; i++) {
                delete[] dataset.data[i];
            }
            delete[] dataset.data;
        }

        if (harmonic_output.height > 0) {
            for (int i = 0; i < harmonic_output.height; i++) {
                delete[] harmonic_output.data[i];
            }
            delete[] harmonic_output.data;
        }
        fftw_free(freq_spect);

    };

    inline void write_harmonics_to_csv(std::string file_dir) const {
        MyFourierClass::write_to_csv<double>(file_dir, this->harmonic_output.data, (int)this->harmonic_output.width, (int)this->harmonic_output.height);
    };

    //Template functions have to be inline (i think)
    template<class T>
    static inline void write_to_csv(std::string file_dir, T** arr, const int width, const int height) {
        std::ofstream out(file_dir + ".csv");

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
    static inline void write_to_csv_1d(const char* file_dir, const T* const arr, const int L) {
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

private:
    myMatrix<double> dataset; // Stores the matrix to be Fourier transformed.
    double Fs; //  sample rate 
    myMatrix<double> harmonic_output; // each row is one harmonic. 
    fftw_complex* freq_spect; // keep this representation for inverse_fft.

    // contains lists of wave parameters. use myhelpers::maxk on amplitudeList to get the most important ones. 
    std::vector<double> amplitudeList;
    std::vector<double> phaseList;
    std::vector<std::complex<double>> frequencyList; // Same as freq_spect but stored less stupidly 
    

    static void forward_fft(const int bins, const size_t L, myMatrix<double> dataset, fftw_complex* out);
    inline void execute_forward_fft(int bins);

    static void fourier_series(const fftw_complex* freq_spect, const int terms, const double Fs, const size_t L, double** out_sin, double** out_cos);

    static void synthesise_from_waves(const int terms, const size_t L,const double* const * sin_mat, const double*const* cos_mat, double** out_synthesis);
    inline void execute_synthesise_from_waves(const int terms, const double* const* out_sin, const double* const* out_cos) {
        this->harmonic_output.data = new double* [terms];
        this->harmonic_output.height = terms;
        this->harmonic_output.width = this->dataset.height;
        synthesise_from_waves(terms, this->dataset.height, out_sin, out_cos, this->harmonic_output.data);
    };

    void inverse_fft(const int terms, fftw_complex* input, double* output) const;
   


};
