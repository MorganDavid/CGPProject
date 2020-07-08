#include "cgp.h"
#include <iostream>
#include <fstream>
#include "mycgp.h"

void cgpWrapper::my_runCGP() {
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;

    int numInputs = 1;
    int numNodes = 15;
    int numOutputs = 1;
    int nodeArity = 2;

    int numGens = 50000;
    int updateFrequency = 500;
    double targetFitness = 0.1;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,mul,sin,cos");

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    printParameters(params);

    trainingData = initialiseDataSetFromFile("complex-300pnts.csv"); // function in https://www.desmos.com/calculator/zlxnnoggsu
    int dssize = getNumDataSetSamples(trainingData);


    chromo = runCGP(params, trainingData, numGens);

    std::ofstream out("predictions.csv");
    for (int i = 0; i < dssize; i++) {
        const double one_input = getDataSetSampleInput(trainingData, i, 0);
        const double arr[] = { one_input };
        executeChromosome(chromo, arr);
        out << one_input << "," << (getChromosomeOutput(chromo, 0)) << "," << getDataSetSampleOutput(trainingData, i, 0) << "," << std::endl;

    }
    out.close();


    system("gnuplot -e \"set datafile separator comma; \
            set terminal jpeg; \
            plot \\\"predictions.csv\\\" using ($1):($2) title \\\"Prediction\\\" with lines,\
            \\\"predictions.csv\\\" using ($1):($3) title \\\"True\\\" with lines;\" > \"results.jpeg\"\"");

    printChromosome(chromo, 0);

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);
}
