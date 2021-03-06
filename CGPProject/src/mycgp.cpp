#include <iostream>
#include <fstream>
#include "mycgp.h"
#include "my_fourier.h"
#include "my_struct_definitions.c"
#include "myhelpers.h"

struct parameters* cgpWrapper::params=NULL;

/*
 > Implementation of approaches are at the bottom of this file

----- Helpers funtions and input functions for approaches ------

*/

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


/*

    Hardcoded hyperparameters for harmonicCGP and CGP.

*/
void cgpWrapper::initializeParams() {
    int numInputs = 1;
    int numNodes = 30;
    int numOutputs = 1;
    int nodeArity = 2;
    const int harmonics_count = 4;

    int updateFrequency = 20; // must be integer multiple of myNumGens
    double targetFitness = 0.1;
    int fourier_terms = 10;
    double** out_synth = new double* [fourier_terms];

    params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity);

    addNodeFunction(params, "add,mul,sin,cos,pi,rand");

    setCustomFitnessFunction(params, MSE, "MSE");

    setTargetFitness(params, targetFitness);
    //setLambda(params, 128);
    //setMu(params, 18);
    setUpdateFrequency(params, updateFrequency);
    //setMutationType(params, "point");
    setMutationRate(params, 0.1);
    setNumThreads(params, 4);
    printParameters(params);
    params->myNumGens = 10000;
    params->myNumRepeats = 1;

    setHarmonicRunParamaters(params, harmonics_count, 0, nullptr);
    setHarmonicRunResultsInit(params, harmonics_count, 50000, updateFrequency);
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
void cgpWrapper::replaceCGPdataSetCol(dataSet* trainingData, const std::vector<double> x, const int col, const bool replaceInput)
{
    for (int i = 0; i < trainingData->numSamples; i++) {
        if (!replaceInput) { trainingData->outputData[i][col] = x[i]; };
        if (replaceInput) { trainingData->inputData[i][col] = x[i]; };
    }
}

// Used to subtract the previous periods best chromo predictions from the next periods target values for the wave like system.
// <returns> target vector with subtraction applied - data not modified. 
std::vector<double> cgpWrapper::subtractChromoPredsFromTarget(struct chromosome* chromo, struct dataSet* data, std::vector<double> new_target) {
    /* for each sample in data */
    for (int i = 0; i < getNumDataSetSamples(data); i++) {
        /* calculate the chromosome outputs for the set of inputs  */
        executeChromosome(chromo, getDataSetSampleInputs(data, i));

        new_target[i] = new_target[i] - getChromosomeOutput(chromo, 0);        
    }
    return new_target;
}

// Subtracts ALL the chromo outputs from data on new target.
// unlike the other version of this function, this function also modifies the `dataSet* data` to represent new_target.
std::vector<double> cgpWrapper::subtractChromoPredsFromTarget(struct chromosome** chromos, const int num_chromos, struct dataSet* data, std::vector<double> new_target) {
    for (int x = 0; x < num_chromos; x++) {
        new_target = subtractChromoPredsFromTarget(chromos[x], data, new_target);
        replaceCGPdataSetCol(data, new_target, 0, false);
    }
    return new_target;
}

void cgpWrapper::writeAndPlot(struct chromosome* chromo, struct dataSet* data, std::string filename) {
    std::string old_filename = filename;
    filename = "output/"+filename+".csv";
    std::ofstream out(filename);
    for (int i = 0; i < data->numSamples; i++) {
        const double one_input = getDataSetSampleInput(data, i, 0);
        const double arr[] = { one_input };
        executeChromosome(chromo, getDataSetSampleInputs(data,i));
        out << one_input << "," << (getChromosomeOutput(chromo, 0)) << "," << getDataSetSampleOutput(data, i, 0) << "," << std::endl;
    }

    out.close();
    std::string one   ("gnuplot -e \"set datafile separator comma; set terminal jpeg; plot \\\"");
    std::string two   ("\\\" using ($1):($2) title \\\"Prediction\\\" with lines,\\\"");
    std::string three ("\\\" using ($1):($3) title \\\"True\\\" with lines;\" > \"");
    std::string four  (old_filename+ ".jpeg\\\"\"");
    std::string cmd   (one+filename+two+filename+three+four);
    system            (cmd.c_str());
    std::cout << cmd << std::endl;
}
// If a list of chromosomes is passed, add up all the chromosomes. used for wave-like aproach. 
void cgpWrapper::writeAndPlot(struct chromosome** chromos, const int num_of_chromos, struct dataSet* data, std::string filename) {
    std::string old_filename = filename;
    filename = "output/" + filename + ".csv";
    std::ofstream out(filename);
    for (int i = 0; i < data->numSamples; i++) {
        const double one_input = getDataSetSampleInput(data, i, 0);
        const double arr[] = { one_input };

        double this_output = 0;
        for (int x = 0; x < num_of_chromos; x++) {
            executeChromosome(chromos[x], getDataSetSampleInputs(data, i));
            this_output += getChromosomeOutput(chromos[x], 0);
        }
        out << one_input << "," << (this_output) << "," << getDataSetSampleOutput(data, i, 0) << "," << std::endl;
    }

    out.close();
    std::string one("gnuplot -e \"set datafile separator comma; set terminal jpeg; plot \\\"");
    std::string two("\\\" using ($1):($2) title \\\"Prediction\\\" with lines,\\\"");
    std::string three("\\\" using ($1):($3) title \\\"True\\\" with lines;\" > \"");
    std::string four(old_filename + ".jpeg\\\"\"");
    std::string cmd(one + filename + two + filename + three + four);
    system(cmd.c_str());
    std::cout << cmd << std::endl;
}

/// <summary>
/// Adds the wave properties of the top k harmonics to the input of `data`. 
/// </summary>
/// <param name="k"> number of harmoncis to add into the CGP inputs. </param>
/// <param name="data"></param>
/// <param name="amplitudeList"></param>
/// <param name="frequencyList"></param>
/// <param name="outData"></param>
struct dataSet* cgpWrapper::constructDataSetWithWaveProperties( const int k, const int Fs, struct dataSet* data, std::vector<double> amplitudeList, std::vector<std::complex<double>> frequencyList ) {
    if (k <= 0) return data;
    
    //get maxk amplitudes
    std::vector<int> maxInds = myhelpers::maxk(amplitudeList, k);
    
    double** inputs = new double*[data->numSamples];
    int newNumInputs = data->numInputs + k * 2; // *3 because we are storing 3 extra inputs into each sample or 3 harmoncis 
    
    for (int i = 0; i < data->numSamples; i++) {
        // Space for the two new inputs.
        inputs[i] = new double[newNumInputs];
        // Fill existing inputs into outData
        for (int x = 0; x < data->numInputs; x++) {
            inputs[i][x] = data->inputData[i][x];
        }

        // Add two new inputs
        for (int y = 0; y < k; y++) {
            inputs[i][data->numInputs + y] = amplitudeList[maxInds[y]];
            inputs[i][data->numInputs + y + k] = (Fs*maxInds[y])/data->numSamples;
        }

    }
    //print inputs to csv for debug 
    // self-init outData from malloc array. 

    //MyFourierClass::write_to_csv("inputs", t->inputData, newNumInputs, data->numSamples);
    //saveDataSet(outData, "ds_with_Fourier_info.csv");
    return initialiseDataSetFromMatrix(newNumInputs, data->numOutputs, data->numSamples, inputs, data->outputData);
}   


/// <summary>
/// Adds the top k harmonic synthesis outputs to the input.
/// It adds the entire harmonic synthesis and not just wave properties
/// </summary>
/// <param name="k"> number of harmoncis to add into the CGP inputs. </param>
/// <param name="data"></param>
/// <param name="amplitudeList"></param>
/// <param name="frequencyList"></param>
/// <param name="outData"></param>
struct dataSet* cgpWrapper::constructDataSetWithFourierInputs(const int k, const int Fs, struct dataSet* data, const double * const * cos_mat, const double * const * sin_mat) {
    if (k <= 0) return data;
    
    double** inputs = new double* [data->numSamples];
    int newNumInputs = data->numInputs + k * 2;

    for (int i = 0; i < data->numSamples; i++) {
        // Space for the two new inputs.
        inputs[i] = new double[newNumInputs];
        // Fill existing inputs into inputData
        for (int x = 0; x < data->numInputs; x++) {
            inputs[i][x] = data->inputData[i][x];
        }

        // Add k new inputs
        for (int y = 0; y < k; y++) {
            inputs[i][data->numInputs + y] = cos_mat[y][i];
            inputs[i][data->numInputs + y + k] = sin_mat[y][i];
        }

    }
    //print inputs to csv for debug 
    // self-init outData from malloc array. 

    MyFourierClass::write_to_csv("inputs_with_fourier.csv", inputs, newNumInputs, data->numSamples);
   // saveDataSet(outData, "ds_with_Fourier_info.csv");
    return initialiseDataSetFromMatrix(newNumInputs, data->numOutputs, data->numSamples, inputs, data->outputData);
}

/*


----- Functions for running the different approaches ------


*/

// Standard CGP (no epochs) but with wave properties data. 
void cgpWrapper::harmonic_runCGP_with_fourier_input(std::string filename) {
    const int harmonics_count = params->harmonicRunParamters->numPeriods;

    struct dataSet* original_data = initialiseDataSetFromFile(filename.c_str());
    struct dataSet* trainingData = initialiseDataSetFromFile(filename.c_str());
    int Fs = 1;
    MyFourierClass f(Fs, original_data);
    f.execute_extract_harmonics(harmonics_count);

    // Modifying data to contain wave properties and input data
    // trainingData will change with each harmonic update, whereas orinal data will always be original input function
    // in input 0.
    original_data = constructDataSetWithWaveProperties(harmonics_count, Fs, original_data, f.get_amplitude_list(), f.get_frequency_list());

    trainingData = constructDataSetWithWaveProperties(harmonics_count, Fs, original_data, f.get_amplitude_list(), f.get_frequency_list());
    setNumInputs(params, trainingData->numInputs);
    saveDataSet(original_data, "original_data.csv");

    // Finally, run on original data
    setHarmonicRunParamaters(params, harmonics_count, harmonics_count - 1, original_data);
    struct chromosome* chromo = my_runCGP(original_data);
    writeAndPlot(chromo, original_data, "plot_final");

    freeDataSet(trainingData);
    freeChromosome(chromo);
    freeParameters(params);

}

// Wave-like system.
void cgpWrapper::harmonic_runCGP_wave(std::string filename) {
    const int harmonics_count = params->harmonicRunParamters->numPeriods;

    struct dataSet* original_data = initialiseDataSetFromFile(filename.c_str());
    struct dataSet* trainingData = initialiseDataSetFromFile(filename.c_str());
    int Fs = 100;
    MyFourierClass f(Fs, original_data);

    // TODO: terms_arr needs to be able to be passed from outside mycgp.
    int* terms_arr = new int[harmonics_count] {1, 2, 4, 8};
    f.set_terms_array(harmonics_count, terms_arr);

    f.execute_extract_harmonics(harmonics_count);

    f.write_harmonics_to_csv("harmonics.csv");
    f.write_to_csv_1d("amps.csv ", &f.get_amplitude_list()[0], f.get_amplitude_list().size());
    f.write_to_csv_1d("freq.csv", &f.get_frequency_list()[0], f.get_frequency_list().size());

    // Modifying data to contain wave properties and input data
    // trainingData will change with each harmonic update, whereas orinal data will always be original input function
    // in input 0.
    original_data = constructDataSetWithWaveProperties(harmonics_count, Fs, original_data, f.get_amplitude_list(), f.get_frequency_list());
    trainingData = constructDataSetWithWaveProperties(harmonics_count, Fs, trainingData, f.get_amplitude_list(), f.get_frequency_list());

    setNumInputs(params, trainingData->numInputs);
    struct chromosome** best_chromos = new struct chromosome* [harmonics_count](); // () initilizes to 0

    // Runs on all harmnonics except the last one. 
    for (int i = 1; i < harmonics_count; i++) {
        // Update trainingData with harmonic data instead.  
        std::vector<double> x = f.getSynthesisWithHarmonics(i);
        if (i != 1) { // Subtract the previoious period best chromo predicitions from this period's target value.
            x = subtractChromoPredsFromTarget(best_chromos, i - 1, trainingData, x);
        }
        else {// if it's the first run, just run on first harmonic.
            replaceCGPdataSetCol(trainingData, x, 0, false);
        }

        saveDataSet(trainingData, std::string("trainingData" + std::to_string(i) + ".csv").c_str());
        //Run CGP on the updated dataset
        setHarmonicRunParamaters(params, harmonics_count, i - 1, original_data);
        best_chromos[i - 1] = my_runCGP(trainingData);

        writeAndPlot(best_chromos, i, trainingData, "plot_" + std::to_string(i));
        std::cout << "Finished Harmonic " << std::endl;
        saveChromosome(best_chromos[i - 1], (std::string("output/chromo") + std::to_string(i)).c_str());
    }

    auto x = subtractChromoPredsFromTarget(best_chromos[harmonics_count - 2], trainingData, myhelpers::getColFromMatrix(original_data->outputData, original_data->numSamples, 0));
    replaceCGPdataSetCol(trainingData, x, 0, false);

    // Finally, run on original data
    setHarmonicRunParamaters(params, harmonics_count, harmonics_count - 1, original_data);
    best_chromos[harmonics_count - 2] = my_runCGP(original_data);
    writeAndPlot(best_chromos[harmonics_count - 2], original_data, "plot_final");

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
// Inter-epoch population pres approach. 
// Change input approach by changing hwo you run original_data and trainingData. 
// Removing line `setInitChromo(params, best_chromos[i - 1]);` prevents parents being carried to next epoch.
void cgpWrapper::harmonic_runCGP(std::string filename) {
    const int harmonics_count = params->harmonicRunParamters->numPeriods;

    struct dataSet* original_data = initialiseDataSetFromFile(filename.c_str()); // function in https://www.desmos.com/calculator/zlxnnoggsu
    struct dataSet* trainingData = initialiseDataSetFromFile(filename.c_str());
    int Fs = 100;
    MyFourierClass f(Fs, original_data);
    f.execute_extract_harmonics(harmonics_count);

    f.write_harmonics_to_csv("harmonics.csv");
    f.write_to_csv_1d("amps.csv ", &f.get_amplitude_list()[0], f.get_amplitude_list().size());
    f.write_to_csv_1d("freq.csv", &f.get_frequency_list()[0], f.get_frequency_list().size());

    // Modifying data to contain wave properties and input data
    // trainingData will change with each harmonic update, whereas orinal data will always be original input function
    // in input 0.
    //original_data = constructDataSetWithWaveProperties(harmonics_count, Fs, original_data, f.get_amplitude_list(), f.get_frequency_list());
    //trainingData = constructDataSetWithWaveProperties(harmonics_count, Fs, trainingData, f.get_amplitude_list(), f.get_frequency_list());

    setNumInputs(params, trainingData->numInputs);
    struct chromosome** best_chromos = new struct chromosome* [harmonics_count]();// () initilizes to 0    

    //Runs on all harmnonics except the last one. 
    for (int i = 1; i < harmonics_count; i++) {
        // Update trainingData with harmonic data instead.  
        std::vector<double> x = f.getSynthesisWithHarmonics(i);

        replaceCGPdataSetCol(trainingData, x, 0, false);
        //Run CGP on the updated dataset
        setHarmonicRunParamaters(params, harmonics_count, i - 1, original_data);
        best_chromos[i - 1] = my_runCGP(trainingData);
        setInitChromo(params, best_chromos[i - 1]);
        writeAndPlot(best_chromos[i - 1], trainingData, "plot_" + std::to_string(i));
        std::cout << "Finished Harmonic " << std::endl;
        saveChromosome(best_chromos[i - 1], (std::string("output/chromo") + std::to_string(i)).c_str());
    }
    // Finally, run on original data
    setHarmonicRunParamaters(params, harmonics_count, harmonics_count - 1, original_data);
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
