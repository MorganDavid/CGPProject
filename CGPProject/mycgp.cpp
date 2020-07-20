#include <iostream>
#include <fstream>
#include "mycgp.h"
#include "my_fourier.h"

void cgpWrapper::harmonic_runCGP() {
    struct parameters* params = NULL;

    int numInputs = 1;
    int numNodes = 20;
    int numOutputs = 1;
    int nodeArity = 2;
    const int harmonics_count = 5;

    int numGens = 5000;
    int updateFrequency = 20;
    double targetFitness = 0.1;
    int fourier_terms = 5;
    double** out_synth = new double* [fourier_terms];

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "rand,add,mul,pi,sin,cos");

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    printParameters(params);

    setHarmonicRunResultsInit(params, harmonics_count, numGens, updateFrequency);

    struct dataSet* original_data = initialiseDataSetFromFile("complex-300pnts.csv"); // function in https://www.desmos.com/calculator/zlxnnoggsu
    struct dataSet* trainingData = initialiseDataSetFromFile("complex-300pnts.csv");

    int dssize = getNumDataSetSamples(trainingData);

    MyFourierClass f(dssize, trainingData);
    f.execute_extract_harmonics(harmonics_count);

    struct chromosome* best_chromos[harmonics_count] = { NULL };

    for (int i = 1; i < 2; i++) {
        // Update trainingData with harmonic data instead.  
        std::vector<double> x = f.getSynthesisWithHarmonics(i);
        replaceCGPdataSetCol(trainingData, x, 0);

        //Run CGP on the updated dataset
        setHarmonicRunParamaters(params, harmonics_count, i - 1, original_data);
        best_chromos[i - 1] = runCGP(params, trainingData, numGens);
        setInitChromo(params, best_chromos[i - 1]);
        std::cout << "Finished Harmonic " << std::endl;
    }
    MyFourierClass::write_to_csv<double>("realfitness", getRealFitnessFromParams(params), numGens / updateFrequency, harmonics_count);
    MyFourierClass::write_to_csv<double>("harmonicfitness", getHarmonicFitnessFromParams(params), numGens / updateFrequency, harmonics_count);

    std::ofstream out("predictions.csv");
    for (int i = 0; i < dssize; i++) {
        const double one_input = getDataSetSampleInput(trainingData, i, 0);
        const double arr[] = { one_input };
        executeChromosome(best_chromos[0], arr);
        out << one_input << "," << (getChromosomeOutput(best_chromos[0], 0)) << "," << getDataSetSampleOutput(trainingData, i, 0) << "," << std::endl;

    }

    out.close();

    system("gnuplot -e \"set datafile separator comma; \
            set terminal jpeg; \
            plot \\\"predictions.csv\\\" using ($1):($2) title \\\"Prediction\\\" with lines,\
            \\\"predictions.csv\\\" using ($1):($3) title \\\"True\\\" with lines;\" > \"results.jpeg\"\"");

    printChromosome(best_chromos[0], 0);

    freeDataSet(trainingData);
    for (int i = 0; i < harmonics_count; i++) {
        if (best_chromos[i] != NULL) {
            freeChromosome(best_chromos[i]);
        }
    }
    freeParameters(params);
}

struct chromosome* cgpWrapper::my_runCGP() {
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;

    int numInputs = 1;
    int numNodes = 10;
    int numOutputs = 1;
    int nodeArity = 2;

    int numGens = 10000;
    int updateFrequency = 500;
    double targetFitness = 0.1;
    int numRepeats = 2;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,pi,mul,sin,1");

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);
    setMutationRate(params, 0.2);
    printParameters(params);
    setNumThreads(params, 4);
    trainingData = initialiseDataSetFromFile("complex-300pnts.csv"); // function in https://www.desmos.com/calculator/zlxnnoggsu
    int dssize = getNumDataSetSamples(trainingData);

    MyFourierClass f(dssize, trainingData);
    f.execute_extract_harmonics(5);

    std::vector<double> x = f.getSynthesisWithHarmonics(1);
    replaceCGPdataSetCol(trainingData, x, 0);
    std::ofstream myfile("trainingData.csv");
    for (int n = 0; n < x.size(); n++)
    {
        myfile << trainingData->outputData[n][0] << std::endl;
    }

    struct results* rels;
    
    //Do cgp
    rels = repeatCGP(params, trainingData, numGens, numRepeats);
    struct chromosome* best_chromo = getBestChromosomeFromResults(rels);

    std::ofstream out("predictions.csv");
    for (int i = 0; i < dssize; i++) {
        const double one_input = getDataSetSampleInput(trainingData, i, 0);
        const double arr[] = { one_input };
        executeChromosome(best_chromo, arr);
        out << one_input << "," << (getChromosomeOutput(best_chromo, 0)) << "," << getDataSetSampleOutput(trainingData, i, 0) << "," << std::endl;
    }

    out.close();


    system("gnuplot -e \"set datafile separator comma; \
            set terminal jpeg; \
            plot \\\"predictions.csv\\\" using ($1):($2) title \\\"Prediction\\\" with lines,\
            \\\"predictions.csv\\\" using ($1):($3) title \\\"True\\\" with lines;\" > \"results.jpeg\"\"");

    printChromosome(best_chromo, 0);

    freeDataSet(trainingData);
    freeResults(rels);
    freeParameters(params);
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

