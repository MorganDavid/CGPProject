#include "fftw3.h"
#include "fstream"
#include <complex>
#include <vector>
#include "mymatrix.h"

/*
 NOTE: All data which goes into this class must have EVEN LENGTH in the first dimension.
 Class object must be destroyed and reinstantiated every time you want to run on new data otherweise it will leak memory.
*/

class MyFourierClass {
public:
    void execute_extract_harmonics(int terms);
    std::vector<double> getSynthesisWithHarmonics(int num_harmonics);

    void setDatasetFromCGP(struct dataSet* data);
    void load_from_csv(const std::string file_dir);

    MyFourierClass(double _Fs, dataSet* _dataset) :Fs{ _Fs } {
        this->setDatasetFromCGP(_dataset);
        // determine Fs based on difference between first two sample inputs
        this->freq_spect = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dataset.height);
    };

    MyFourierClass(double _Fs, std::string file_dir) {
        this->load_from_csv(file_dir);
        // determine Fs based on difference between first two sample inputs
        this->freq_spect = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->dataset.height);
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
        MyFourierClass::write_to_csv<double>(file_dir + ".csv", this->harmonic_output.data, this->harmonic_output.width, this->harmonic_output.height);
        for (int i = 0; i < (int)this->harmonic_output.height; i++) {
            MyFourierClass::write_to_csv_1d<double>(file_dir+"_"+std::to_string(i)+".csv", this->harmonic_output.data[i], (int)this->harmonic_output.width);
        }
     };

    //Template functions have to be inline (i think)
    template<class T>
    static inline void write_to_csv(std::string file_dir, T** arr, const int width, const int height) {
        std::ofstream out("output/"+file_dir);

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
    //TODO: see if there's a less stupid way to have all these similar functions in one function. 
    template<class T>
    static inline void write_to_csv_1d(const std::string file_dir, const T* const arr, const int L) {
        std::ofstream out("output/"+file_dir);

        for (int r = 0; r < L; r++) {
            out << arr[r] << ',' << '\n';
        }
        out.close();
    };
   
    //TODO: make a new implementation of this to generalize to any type with tostring ability. 
    static inline void write_to_csv(const std::string file_dir, fftw_complex* arr, const int L) {
        std::ofstream out("output/"+file_dir);

        for (int r = 0; r < L; r++) {
            out << arr[r][0] << ":" << arr[r][1] << ',' << "\n";
        }
        out.close();
    };
    inline std::vector<double> get_amplitude_list() {
        return amplitudeList;
    }
    inline std::vector < std::complex<double> > get_frequency_list() {
        return frequencyList;
    }
    inline std::vector<double> get_phase_list() {
        return phaseList;
    }
    /// If we don't want to increment by one for every term, use this. define as [2,6,10,20] for example.
    /// Run this before `execute_extract_harmonics` (optinal)
    inline void set_terms_array(int terms, int* terms_arr) {
        this->num_of_terms = terms;
        this->terms_array = terms_arr;
    }
    inline int* get_terms_array() {
        return terms_array;
    }
private:
    myMatrix<double> dataset; // Stores the matrix to be Fourier transformed.
    int num_of_terms = -1; // Number of terms to generatie 
    int* terms_array = NULL; // If we don't want to increment by one for every term, use this. define as [2,6,10,20] with four num_of_terms for example.
    double Fs; //  sample rate 
    myMatrix<double> harmonic_output;
        
    fftw_complex* freq_spect;

    std::vector<double> amplitudeList;
    std::vector<double> phaseList;
    std::vector<std::complex<double>> frequencyList;

    static void forward_fft(const int bins, const size_t L, myMatrix<double> dataset, fftw_complex* out);
    inline void execute_forward_fft(int bins) {
        forward_fft(bins, this->dataset.height, this->dataset, this->freq_spect);
    };
    static std::vector<double> calculate_amplitude_list(std::vector<std::complex<double>> freq_spect);

	static void fourier_series(const std::vector<std::complex<double>> freq_spect_cmplx, const std::vector<double> amplitudeList, const int terms, const double Fs, const size_t L, std::complex<double>* out_synthesis);


    void inverse_fft(const int terms, fftw_complex* input, double* output) const;


};
