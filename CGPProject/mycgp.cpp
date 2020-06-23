#include "cgp.h"
#include <iostream>

class cgpWrapper {
public:
    static void my_runCGP();
private:
    static double fitness(struct parameters*, struct chromosome*, struct dataSet*);
};

double cgpWrapper::fitness(struct parameters* params, struct chromosome* chromo, struct dataSet* data) {
    int i, j;
    double error = 0;

    /* error checking */
    if (getNumChromosomeInputs(chromo) != getNumDataSetInputs(data)) {
        printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
        printf("Terminating CGP-Library.\n");
        exit(0);
    }

    if (getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)) {
        printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
        printf("Terminating CGP-Library.\n");
        exit(0);
    }

    /* for each sample in data */
    for (i = 0; i < getNumDataSetSamples(data); i++) {

        /* calculate the chromosome outputs for the set of inputs  */
        executeChromosome(chromo, getDataSetSampleInputs(data, i));

        /* for each chromosome output */
        for (j = 0; j < getNumChromosomeOutputs(chromo); j++) {

            error += pow(getChromosomeOutput(chromo, j) - getDataSetSampleOutput(data, i, j), 2);
        }
    }

    return error / getNumDataSetSamples(data);
}

void cgpWrapper::my_runCGP() {
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;

    int numInputs = 1;
    int numNodes = 20;
    int numOutputs = 1;
    int nodeArity = 2;

    int numGens = 100000;
    int updateFrequency = 500;
    double targetFitness = 0.0;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,sin,mult,pi");

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    setCustomFitnessFunction(params, cgpWrapper::fitness, "MSE");

    printParameters(params);

    trainingData = initialiseDataSetFromFile("sinwaves.csv"); // function in sinewaves.csv https://www.desmos.com/calculator/uoafllrqo8
    int dssize = getNumDataSetSamples(trainingData);


    chromo = runCGP(params, trainingData, numGens);

    std::cout << "chromo outputs: " << std::endl;
    std::cout << "x : ypred : ytrue" << std::endl;
    for (int i = 0; i < dssize; i++) {
        const double one_input = getDataSetSampleInput(trainingData, i, 0);
        const double arr[] = { one_input };
        executeChromosome(chromo, arr);
        std::cout << one_input << " : " << (getChromosomeOutput(chromo, 0)) << " : " << getDataSetSampleOutput(trainingData, i, 0) << std::endl;
    }


    printChromosome(chromo, 0);

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);
}
