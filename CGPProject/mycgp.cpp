#include <iostream>
#include <fstream>
#include "mycgp.h"
#include "my_fourier.h"
#include "my_struct_definitions.c"

struct parameters* cgpWrapper::params=NULL;

void cgpWrapper::initializeParams() {
    int numInputs = 1;
    int numNodes = 20;
    int numOutputs = 1;
    int nodeArity = 2;
    const int harmonics_count = 3;

    int updateFrequency = 20; // must be integer multiple of myNumGens
    double targetFitness = 0.1;
    int fourier_terms = 3;
    double** out_synth = new double* [fourier_terms];

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "1,div,add,mul,pi,sin,cos");

    setTargetFitness(params, targetFitness);
    setLambda(params, 5);
    setUpdateFrequency(params, updateFrequency);
    setNumThreads(params, 4);
    printParameters(params);
    params->myNumGens = 5000;
    params->myNumRepeats = 5;

    setHarmonicRunParamaters(params, harmonics_count, 0, nullptr);
    setHarmonicRunResultsInit(params, harmonics_count, 30000, updateFrequency);
}

void cgpWrapper::harmonic_runCGP() {
    const int harmonics_count = params->harmonicRunParamters->numPeriods;

    struct dataSet* original_data = initialiseDataSetFromFile("complex-300pnts.csv"); // function in https://www.desmos.com/calculator/zlxnnoggsu
    struct dataSet* trainingData = initialiseDataSetFromFile("complex-300pnts.csv");

    int dssize = getNumDataSetSamples(trainingData);

    MyFourierClass f(dssize, trainingData);
    f.execute_extract_harmonics(harmonics_count);

    struct chromosome** best_chromos = new struct chromosome*[harmonics_count]();// () initilizes to 0    

    for (int i = 1; i <= harmonics_count; i++) {
        // Update trainingData with harmonic data instead.  
         std::vector<double> x = f.getSynthesisWithHarmonics(i);
         replaceCGPdataSetCol(trainingData, x, 0);
         if (i == 3) params->myNumGens = 30000;
        //Run CGP on the updated dataset
        setHarmonicRunParamaters(params, harmonics_count, i - 1, original_data);
        best_chromos[i - 1] = my_runCGP(trainingData);
        setInitChromo(params, best_chromos[i - 1]);
        writeAndPlot(best_chromos[i - 1], trainingData, "plot_" + std::to_string(i));
        std::cout << "Finished Harmonic " << std::endl;
    }

    MyFourierClass::write_to_csv<double>("realfitness", getRealFitnessFromParams(params), params->myNumGens / params->updateFrequency, harmonics_count);
    MyFourierClass::write_to_csv<double>("harmonicfitness", getHarmonicFitnessFromParams(params), params->myNumGens / params->updateFrequency, harmonics_count);

    freeDataSet(trainingData);
    for (int i = 0; i < harmonics_count; i++) {
        if (best_chromos[i] != 0) {
            freeChromosome(best_chromos[i]);
        }
    }
    freeParameters(params);
}

struct chromosome* cgpWrapper::my_runCGP(struct dataSet* trainingData) {
    //Do cgp
    struct chromosome* best_chromo = getBestChromosomeFromResults(repeatCGP(params, trainingData, params->myNumGens, params->myNumRepeats));

    printChromosome(best_chromo, 0);

    return best_chromo;
}

/// <summary>
/// Replaces col in the outputData of the dataset with the data stored in the vector x. O(n) operation! 
/// </summary>
/// <param name="trainingData"></param>
/// <param name="x"></param>
/// <param name="col"></param>
void cgpWrapper::replaceCGPdataSetCol(dataSet* trainingData, const std::vector<double> x, const int col)
{
    for (int i = 0; i < trainingData->numSamples; i++) {
        trainingData->outputData[i][col] = x[i];
    }
}

void cgpWrapper::writeAndPlot(struct chromosome* chromo, struct dataSet* data, std::string filename) {
    std::string old_filename = filename;
    filename = filename+".csv";
    std::ofstream out(filename);
    for (int i = 0; i < data->numSamples; i++) {
        const double one_input = getDataSetSampleInput(data, i, 0);
        const double arr[] = { one_input };
        executeChromosome(chromo, arr);
        out << one_input << "," << (getChromosomeOutput(chromo, 0)) << "," << getDataSetSampleOutput(data, i, 0) << "," << std::endl;
    }

    out.close();
    std::string one   ("gnuplot -e \"set datafile separator comma; set terminal jpeg; plot \\\"");
    std::string two   ("\\\" using ($1):($2) title \\\"Prediction\\\" with lines,\\\"");
    std::string three ("\\\" using ($1):($3) title \\\"True\\\" with lines;\" > \"");
    std::string four  (old_filename+ ".jpeg\\\"\"");
    std::string cmd   (one+filename+two+filename+three+four);
    system            (cmd.c_str());
}