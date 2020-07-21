#include <vector>
#include "cgp_harms.h"

class cgpWrapper {
public:
    static void harmonic_runCGP();
    static struct chromosome* my_runCGP(struct dataSet*);
    static void initializeParams();
    static struct parameters* params;
private:
    static void replaceCGPdataSetCol(dataSet* trainingData,const std::vector<double> x,const int col);
    static void writeAndPlot(chromosome* chromo, dataSet* data, std::string filename);
};