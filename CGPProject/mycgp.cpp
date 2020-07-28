#include <iostream>
#include <fstream>
#include "mycgp.h"
#include "my_fourier.h"
#include "my_struct_definitions.c"
#include "myhelpers.h"

struct parameters* cgpWrapper::params=NULL;

double cgpWrapper::MSE(struct parameters* params, struct chromosome* chromo, struct dataSet* data) {

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

            error += pow(getChromosomeOutput(chromo, j) - getDataSetSampleOutput(data, i, j),2);
        }
    }
    error /= getNumDataSetSamples(data);
    return error;
}

void cgpWrapper::initializeParams() {
    int numInputs = 1;
    int numNodes = 18;
    int numOutputs = 1;
    int nodeArity = 2;
    const int harmonics_count = 3;

    int updateFrequency = 20; // must be integer multiple of myNumGens
    double targetFitness = 0.001;
    int fourier_terms = 3;
    double** out_synth = new double* [fourier_terms];

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,mul,pi,sin,cos,rand");

    setCustomFitnessFunction(params, MSE, "MSE");
    
    setTargetFitness(params, targetFitness);
   // setLambda(params, 4);
   // setMu(params, 1);
    setUpdateFrequency(params, updateFrequency);
    setMutationRate(params, 0.6);
    setNumThreads(params, 4);
    printParameters(params);
    params->myNumGens = 10000;
    params->myNumRepeats = 3;

    setHarmonicRunParamaters(params, harmonics_count, 0, nullptr);
    setHarmonicRunResultsInit(params, harmonics_count, 50000, updateFrequency);
}

void cgpWrapper::harmonic_runCGP(std::string filename) {
    const int harmonics_count = params->harmonicRunParamters->numPeriods;

    struct dataSet* original_data = initialiseDataSetFromFile(filename.c_str()); // function in https://www.desmos.com/calculator/zlxnnoggsu
    struct dataSet* trainingData = initialiseDataSetFromFile(filename.c_str());

    int dssize = getNumDataSetSamples(trainingData);

    MyFourierClass f(dssize, trainingData);
    f.execute_extract_harmonics(harmonics_count);
    struct dataSet* dsWithExtraInputs = nullptr;
    constructDataset(3, original_data, f.get_amplitude_list(), f.get_frequency_list(), dsWithExtraInputs);
    struct chromosome** best_chromos = new struct chromosome*[harmonics_count]();// () initilizes to 0    

    //Runs on all harmnonics except the last one. 
    for (int i = 1; i < harmonics_count; i++) {
        // Update trainingData with harmonic data instead.  
         std::vector<double> x = f.getSynthesisWithHarmonics(i);
         replaceCGPdataSetCol(trainingData, x, 0);
        //Run CGP on the updated dataset
        setHarmonicRunParamaters(params, harmonics_count, i - 1, original_data);
        best_chromos[i - 1] = my_runCGP(trainingData);
        setInitChromo(params, best_chromos[i - 1]);
        writeAndPlot(best_chromos[i - 1], trainingData, "plot_" + std::to_string(i));
        std::cout << "Finished Harmonic " << std::endl;
    }

    // Finally, run on original data
    params->myNumGens = 100000;
    setHarmonicRunParamaters(params, harmonics_count, harmonics_count-1, trainingData);
    best_chromos[harmonics_count - 1] = my_runCGP(original_data);
    writeAndPlot(best_chromos[harmonics_count - 1], original_data, "plot_final");

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
    writeAndPlot(best_chromo, trainingData, "singlerun");
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
/// <summary>
/// 
/// </summary>
/// <param name="k"> number of harmoncis to add into the CGP inputs. </param>
/// <param name="data"></param>
/// <param name="amplitudeList"></param>
/// <param name="frequencyList"></param>
/// <param name="outData"></param>
void cgpWrapper::constructDataset( const int k, struct dataSet* data, std::vector<double> amplitudeList, std::vector<std::complex<double>> frequencyList, struct dataSet* outData ) {
    if (k <= 0) {
        outData = data;
        return;
    }
    
    //get maxk amplitudes
    std::vector<int> maxInds = myhelpers::maxk(amplitudeList, k);
    
    double** inputs = new double*[data->numSamples];
    int newNumInputs = data->numInputs + 1+k*3; // *3 because we are storing 3 extra inputs into each sample. 
    
    for (int i = 0; i < data->numSamples; i++) {
        // Space for the two new inputs.
        inputs[i] = new double[newNumInputs];
        // Fill existing inputs into outData
        for (int x = 0; x < data->numInputs; x++) {
            inputs[i][x] = data->inputData[i][x];
        }
        // Add two new inputs

        for (int y = 0; y < k; y++) {
            inputs[i][data->numInputs + y] = frequencyList[maxInds[y]].imag(); 
            inputs[i][data->numInputs + y + 3] = amplitudeList[maxInds[y]];
            inputs[i][data->numInputs + y + 6] = frequencyList[maxInds[y]].real();
        }

    }

    outData = initialiseDataSetFromArrays(newNumInputs, data->numOutputs, data->numSamples, inputs[0], data->outputData[0]);
    saveDataSet(outData, "ds_with_Fourier_info.csv");
}   