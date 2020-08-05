#include <vector>
#include "cgp_harms.h"
#include <complex>

class cgpWrapper {
public:
    static void harmonic_runCGP(std::string);
    static struct chromosome* my_runCGP(struct dataSet*);
    static void initializeParams();
    static double MSE(struct parameters* params, struct chromosome* chromo, struct dataSet* data);
    static struct parameters* params;
private:
    static void replaceCGPdataSetCol(dataSet* trainingData,const std::vector<double> x,const int col);
    static void writeAndPlot(chromosome* chromo, dataSet* data, std::string filename);
	static  dataSet* constructDataSetWithWaveProperties(const int k, const int Fs, dataSet* data, std::vector<double> amplitudeList, std::vector<std::complex<double>> frequencyList);
};