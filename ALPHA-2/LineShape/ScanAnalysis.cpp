#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/Algorithms.h"
using namespace RooFit;

double mean(std::vector<double> v){
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	return sum/ v.size();
}

double stdev(std::vector<double> v){
	double accum = 0.0;
	double m = mean(v);
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
    	accum += (d - m) * (d - m);
	});
	return sqrt(accum/(v.size() - 1));
}

std::vector<double> ScanAnalysis(TString directory = "linear/",
					int start = 0,
					int stop = 999,
					double N = 1, 			// Threshold coefficient
					double fraction = 0.1 	// constant fraction discrimination
					){
	TString ConfFile = directory + "ToyConfiguration.txt";
	std::cout << ConfFile << std::endl;
	//gInterpreter->GenerateDictionary("ToyParser","../Headers/ConfigurationParser.h");
	//gInterpreter->GenerateDictionary("ReadFiles", "../Headers/AnalysisLineShape.h");
	
	ReadConfFile Params(ConfFile); // Read the values from the configuration file 
	Params.Print();
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
	double FrequencyStep = Params.FrequencyStep;
	double startPdf1 = Params.x_cb_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.x_da_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep d-a
	double SweepStep = Params.SweepStep; // number of bin for each lineshape
	
	std::vector<std::string> FileList;
	FileList = getFiles(start,stop, directory); // file list to be analyzed

	// Define some vectors to store the results
	vector<double> MCtruth;
	vector<double> v1_2017, v2_2017;
	vector<double> v1_rev, v2_rev;
	vector<double> v1_thr, v2_thr;
	vector<double> v1_cfrac, v2_cfrac;
	vector<double> v1_neigh, v2_neigh;
	vector<double> diff_2017, diff_rev, diff_thr, diff_cfrac, diff_neigh;
	int count = 0;
	// IMPLEMENTING THE TOY FOR THE ALGORITHM
	for(int i = 0; i < FileList.size(); i += 2){
		count  += 1; std::cout << "Analizzo DataFrame " << count << "\n" << std::endl;
		
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});		// Load i-th dataset
		auto Spectra1 = frame.Filter("runNumber == 1")
							 .Filter("frequence <= 1000")
							 //.Filter("type != 2")
							 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
		
		auto Spectra2 = frame.Filter("runNumber == 1")
							 .Filter("frequence >= 1000")
							 //.Filter("type != 2")
							 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");
		
		auto frame1 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Mean<double>("lineShift");
		auto frame2 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Mean<double>("lineShift");
		// Load the shifts of the lineshapes
		auto frame3 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Take<double>("lineShift");
		auto frame4 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Take<double>("lineShift");
		auto lineShiftda = frame3.GetValue(); std::cout << "LineShift c to b : " << lineShiftda[0] << std::endl;
		auto lineShiftcb = frame4.GetValue(); std::cout << "LineShift c to b : " << lineShiftcb[0] << std::endl;
		
		double onset1;				// Reconstructed onset
		double onset2;				// Reconstructed onset
		double threshold = N*CosmicBackground; // threshold considering the cosmic background
		
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back((Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// THRESHOLD
		onset1 = firstOverThreshold(Spectra1, N);
		onset2 = firstOverThreshold(Spectra2, N);
		v1_thr.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_thr.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_thr.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// 2017 FOWARD
		onset1 = algorithm_2017(Spectra1);
		onset2 = algorithm_2017(Spectra2);
		v1_2017.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_2017.push_back(onset2 - (Params.x_da_start + lineShiftda[0])); 
		diff_2017.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// REVERSED 2017 
		onset1 = reverse_2017(Spectra1);
		onset2 = reverse_2017(Spectra2);
		v1_rev.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_rev.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_rev.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// CONSTANT FRACTION
		onset1 = constFrac(Spectra1, fraction);
		onset2 = constFrac(Spectra2, fraction);
		v1_cfrac.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_cfrac.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_cfrac.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// SUM NEIGHBORS
		onset1 = sumNeighbors(Spectra1, CosmicBackground, N);
		onset2 = sumNeighbors(Spectra2, CosmicBackground, N);
		v1_neigh.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_neigh.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_neigh.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
	}
	// THRESHOLD, FOWARD, REVERSED, CONSTANT FRACTION, SUM NEIGHBORS
	return {mean(diff_thr), stdev(diff_thr),
			mean(diff_2017),stdev(diff_2017),
			mean(diff_rev), stdev(diff_rev),
			mean(diff_cfrac), stdev(diff_cfrac),
			mean(diff_neigh), stdev(diff_neigh)};
}

