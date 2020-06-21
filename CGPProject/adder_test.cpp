/*
#include "cpp.h"
#include <iostream>

class cgpWrapper {
public:
    static void my_runCGP();
private:
    static double fullAdder(struct parameters*, struct chromosome*, struct dataSet*);
};
void cgpWrapper::my_runCGP() {
    struct parameters* params = NULL;
    struct dataSet* trainingData = NULL;
    struct chromosome* chromo = NULL;

    int numInputs = 3;
    int numNodes = 5;
    int numOutputs = 2;
    int nodeArity = 2;

    int numGens = 10000;
    int updateFrequency = 500;
    double targetFitness = 0.0;

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,xor,or");

    setTargetFitness(params, targetFitness);

    setUpdateFrequency(params, updateFrequency);

    setCustomFitnessFunction(params, fullAdder, "fullAdder");

    printParameters(params);

    // full adder input output pairs
    double inputs[8][3] = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
    double outputs[8][2] = { {0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1} };

    trainingData = initialiseDataSetFromArrays(numInputs, numOutputs, 8, inputs[0], outputs[0]);

    chromo = runCGP(params, trainingData, numGens);

    printChromosome(chromo, 0);

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);
}

double cgpWrapper::fullAdder(struct parameters* params, struct chromosome* chromo, struct dataSet* dat) {

    int i;
    double error = 0;

    // full adder truth table
    double inputs[8][3] = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
    double outputs[8][2] = { {0,0},{1,0},{1,0},{0,1},{1,0},{0,1},{0,1},{1,1} };

    //for each line in the truth table
    for (i = 0; i < 8; i++) {

        // calculate the chromosome outputs for each set of inputs
        executeChromosome(chromo, inputs[i]);

        // If the first chromosome outputs differ from the correct outputs increment the error
        if (outputs[i][0] != getChromosomeOutput(chromo, 0)) {
            error++;
        }

        // If the second chromosome outputs differ from the correct outputs increment the error
        if (outputs[i][1] != getChromosomeOutput(chromo, 1)) {
            error++;
        }
    }

    return error;
}

*/