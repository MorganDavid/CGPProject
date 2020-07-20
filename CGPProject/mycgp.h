#include <vector>
#include "cgp_harms.h"
class cgpWrapper {
public:
    static void harmonic_runCGP();
    static struct chromosome* my_runCGP();
private:
    static double fitness(struct parameters*, struct chromosome*, struct dataSet*);
    static void replaceCGPdataSetCol(dataSet* trainingData,const std::vector<double> x,const int col);
};