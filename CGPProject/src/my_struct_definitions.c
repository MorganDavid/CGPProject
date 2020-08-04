//Contains structure defintitions that can be used in C and C++ without including cgp.c. 

#define MUTATIONTYPENAMELENGTH 21
#define FUNCTIONSETSIZE 50
#define FUNCTIONNAMELENGTH 11
#define FITNESSFUNCTIONNAMELENGTH 21
#define SELECTIONSCHEMENAMELENGTH 21
#define REPRODUCTIONSCHEMENAMELENGTH 21
struct parameters {
	int mu;
	int lambda;
	char evolutionaryStrategy;
	double mutationRate;
	double recurrentConnectionProbability;
	double connectionWeightRange;
	int numInputs;
	int numNodes;
	int numOutputs;
	int arity;
	struct functionSet* funcSet;
	double targetFitness;
	int updateFrequency;
	int shortcutConnections;
	void (*mutationType)(struct parameters* params, struct chromosome* chromo);
	char mutationTypeName[MUTATIONTYPENAMELENGTH];
	double (*fitnessFunction)(struct parameters* params, struct chromosome* chromo, struct dataSet* dat);
	char fitnessFunctionName[FITNESSFUNCTIONNAMELENGTH];
	void (*selectionScheme)(struct parameters* params, struct chromosome** parents, struct chromosome** candidateChromos, int numParents, int numCandidateChromos);
	char selectionSchemeName[SELECTIONSCHEMENAMELENGTH];
	void (*reproductionScheme)(struct parameters* params, struct chromosome** parents, struct chromosome** children, int numParents, int numChildren);
	char reproductionSchemeName[REPRODUCTIONSCHEMENAMELENGTH];
	int numThreads;

	int myNumGens;
	int myNumRepeats;
	/* Harmonics stuff */
	//If harmonicRunParameters is nonnull, then we are doing a harmonic run. 

	 // If this is non NULL, initialise the population with copies of this chromo, not random generation.
	// params obj is NOT responsible for initChromo
	struct chromosome* initChromo;
	struct harmonicRunResults* harmonicRunResults; // Stores results when running in harmonic mode. 
	struct harmonicRunParameters* harmonicRunParamters;
};

struct chromosome {
	int numInputs;
	int numOutputs;
	int numNodes;
	int numActiveNodes;
	int arity;
	struct node** nodes;
	int* outputNodes;
	int* activeNodes;
	double fitness;
	double* outputValues;
	struct functionSet* funcSet;
	double* nodeInputsHold;
	int generation;
};

struct dataSet {
	int numSamples;
	int numInputs;
	int numOutputs;
	double** inputData;
	double** outputData;
};
struct node {
	int function;
	int* inputs;
	double* weights;
	int active;
	double output;
	int maxArity;
	int actArity;
};

struct functionSet {
	int numFunctions;
	char functionNames[FUNCTIONSETSIZE][FUNCTIONNAMELENGTH];
	int maxNumInputs[FUNCTIONSETSIZE];
	double (*functions[FUNCTIONSETSIZE])(const int numInputs, const double* inputs, const double* connectionWeights);
};

struct harmonicRunResults {
	double** harmonicFitness;
	double** realFitness;
};

struct harmonicRunParameters {
	int numPeriods;
	int currentPeriod;
	struct dataSet* rawDataset;//Raw dataset. Doesnt' get fourier transformed. Use for real fitness.	
};

struct results {
	int numRuns;
	struct chromosome** bestChromosomes;
};
