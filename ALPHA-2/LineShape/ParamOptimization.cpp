#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/Algorithms.h"
using namespace RooFit;

struct Result {
	vector<double> Mean_Sigma;
	vector<double> Mean_SquareResidual;
} ;

std::vector<double> ParamOptimization(	TString directory,
					TString ConfFile,	// Configuration files
					int stop,		// Number of runs to be analysed
					double Nfilter,		// Filter for the running sum
					double fraction,	// Parameter of Constant Fraction
					double Nsigma,		// Paratemer for t-student test
					double Nthr,		// Parameter of Threshold algorithm
					double thr1,
					double thr2,
					TString folder		// Where to save the plots
					){
	//ConfFile = directory + "ToyConfiguration.txt";
	//std::cout << ConfFile << std::endl;

	ReadConfFile Params(ConfFile); // Read the values from the configuration file 
	Params.Print();
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
	double FrequencyStep = Params.FrequencyStep;
	double startPdf1 = Params.cb_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.da_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep d-a
	
	std::vector<std::string> FileList;
	FileList = getFiles(0,stop, directory); // file list to be analyzed

	// Define some vectors to store the results
	vector<double> MCtruth;
	vector<double> v1bk_2017, v2bk_2017;
	vector<double> v1bk_rev, v2bk_rev;
	vector<double> v1bk_thr, v2bk_thr;
	vector<double> v1_cfrac, v2_cfrac;
	vector<double> v1_neigh, v2_neigh;
	vector<double> v1_sign, v2_sign;
	vector<double> diff_bk_2017;
	vector<double> diff_bk_rev;
	vector<double> diff_bk_thr;
	vector<double> diff_cfrac;
	vector<double> diff_neigh;
	vector<double> diff_sign;
	
	vector<double> x_2017;
	vector<double> x_rev;
	vector<double> x_thr;
	vector<double> x_cfrac;
	vector<double> x_neigh;
	vector<double> x_sign;
	
	// Define the struct for the return
	Result object;
	
	
	int count = 0;
	// IMPLEMENTING THE TOY FOR THE ALGORITHM
	for(int i = 0; i < FileList.size(); i += 2){
		count  += 1; 
		//std::cout << directory << " Analizzo DataFrame " << count << std::endl;
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});		// Load i-th dataset
		
		// C TO B TRANSITION
		// Extract from the dataset the starting frequency of the transition
		auto sweepStart_cb = frame.Filter("repetition == 0")
		                     .Filter("mwfrequence <= 1000")   // filter on the micro wave frequency to identify the c to b transition
		                     .Take<double>("sweepStart");
		auto actualStart_cb = sweepStart_cb.GetValue();       // get starting frequency from Rdataframe node
		
		// Load the data and save it in a Histogram (raw lineshape)
		// the histogram start from actualStart_cb, with a step given by Params.FrequencyStep
		
		auto Spectra1 = frame.Filter("repetition == 0")
			 .Filter("mwfrequence <= 1000")
			 //.Filter("type != 2")
			 .Histo1D({"Counts","Frequence", static_cast<int>(Params.SweepStep),actualStart_cb[0], actualStart_cb[0] + Params.SweepStep*Params.FrequencyStep }, "mwfrequence");
		
		// D TO A TRANSITION
		// Extract from the dataset the starting frequency of the transition	
		auto sweepStart_da = frame.Filter("repetition == 0")
		                     .Filter("mwfrequence >= 1000")
		                     .Take<double>("sweepStart");
		auto actualStart_da = sweepStart_da.GetValue(); // get starting frequency from Rdataframe node
		
		// Load the data and save it in a Histogram (raw lineshape)
		 // the histogram start from actualStart_da, with a step given by Params.FrequencyStep
		auto Spectra2 = frame.Filter("repetition == 0")
			 .Filter("mwfrequence >= 1000")
			 //.Filter("type != 2")
			 .Histo1D({"Counts","Frequence", static_cast<int>(Params.SweepStep), actualStart_da[0], actualStart_da[0]+ Params.SweepStep*Params.FrequencyStep}, "mwfrequence");

		// Load the shifts of the lineshapes
		auto getOnsetCB = frame.Filter("repetition == 0").Filter("mwfrequence >= 1000").Take<double>("trueOnset");
		auto getOnsetDA = frame.Filter("repetition == 0").Filter("mwfrequence <= 1000").Take<double>("trueOnset");
		auto onsetda = getOnsetDA.GetValue();
		auto onsetcb = getOnsetCB.GetValue();
		
		double onset1;				// Reconstructed onset
		double onset2;				// Reconstructed onset
		double threshold = Nthr*CosmicBackground; // threshold considering the cosmic background
		
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back((onsetda[0] - onsetcb[0]));
		
		// THRESHOLD
			// With Cosmic Background
		onset1 = firstOverThreshold(Spectra1, Nthr, CosmicBackground);
		onset2 = firstOverThreshold(Spectra2, Nthr, CosmicBackground);
		v1bk_thr.push_back(onset1 - onsetcb[0]);
		v2bk_thr.push_back(onset2 - onsetda[0]);
		x_thr.push_back(onset2 - onset1);
		diff_bk_thr.push_back(onset2 - onset1 - (onsetda[0] - onsetcb[0]));		
		
		
		// 2017 FOWARD
			// With Cosmic Background
		onset1 = algorithm_2017(Spectra1, CosmicBackground, thr1, thr2);
		onset2 = algorithm_2017(Spectra2, CosmicBackground, thr1, thr2);
		v1bk_2017.push_back(onset1 - onsetcb[0]);
		v2bk_2017.push_back(onset2 - onsetda[0]);
		x_2017.push_back(onset2 - onset1);
		diff_bk_2017.push_back(onset2 - onset1 - (onsetda[0] - onsetcb[0]));

		
		// REVERSED 2017
			// with cosmic background
		onset1 = reverse_2017(Spectra1, CosmicBackground, thr1, thr2);
		onset2 = reverse_2017(Spectra2, CosmicBackground, thr1, thr2);
		v1bk_rev.push_back(onset1 - onsetcb[0]);
		v2bk_rev.push_back(onset2 - onsetda[0]);
		x_rev.push_back(onset2 - onset1);
		diff_bk_rev.push_back(onset2 - onset1 - (onsetda[0] - onsetcb[0]));
		
		
		// CONSTANT FRACTION
			// with cosmic subtraction
		onset1 = constFrac(Spectra1, fraction, CosmicBackground,Nfilter);
		onset2 = constFrac(Spectra2, fraction, CosmicBackground,Nfilter);
		v1_cfrac.push_back(onset1 - onsetcb[0]);
		v2_cfrac.push_back(onset2 - onsetda[0]);
		x_cfrac.push_back(onset2 - onset1);
		diff_cfrac.push_back(onset2 - onset1 - (onsetda[0] - onsetcb[0]));
		
		
		// SUM NEIGHBORS
		onset1 = sumNeighbors(Spectra1, Nsigma, CosmicBackground, Nfilter);
		onset2 = sumNeighbors(Spectra2, Nsigma, CosmicBackground, Nfilter);
		v1_neigh.push_back(onset1 - onsetcb[0]);
		v2_neigh.push_back(onset2 - onsetda[0]);
		x_neigh.push_back(onset2 - onset1);
		diff_neigh.push_back(onset2 - onset1 - (onsetda[0] - onsetcb[0]));
		
		// SIGNIFICANCE
		onset1 = Significance(Spectra1, Nsigma, CosmicBackground, Nfilter);
		onset2 = Significance(Spectra2, Nsigma, CosmicBackground, Nfilter);
		v1_sign.push_back(onset1 - onsetcb[0]);
		v2_sign.push_back(onset2 - onsetda[0]);
		x_sign.push_back(onset2 - onset1);
		diff_sign.push_back(onset2 - onset1 - (onsetda[0] - onsetcb[0]));
		} // Loop on Runs
	
	
	
	return 		{mean(diff_bk_thr),	stdev(diff_bk_thr),		//with background
			mean(diff_bk_2017),	stdev(diff_bk_2017),		//with background
			mean(diff_bk_rev),	stdev(diff_bk_rev),		//with background
			mean(diff_cfrac),	stdev(diff_cfrac),       	//with backgrounf
			mean(diff_neigh),	stdev(diff_neigh),      	// SUM NEIGHBORS
			mean(diff_sign),	stdev(diff_sign),
			stdev(x_thr),	Corr(x_thr,MCtruth),
			stdev(x_2017),	Corr(x_2017,MCtruth),
			stdev(x_rev),	Corr(x_rev,MCtruth),
			stdev(x_cfrac),	Corr(x_cfrac,MCtruth),
			stdev(x_neigh),	Corr(x_neigh,MCtruth),
			stdev(x_sign),	Corr(x_sign, MCtruth),
			stdev(MCtruth)};	
}

