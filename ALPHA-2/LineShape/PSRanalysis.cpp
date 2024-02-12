#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"

//#include <Fit/Fitter.h>
//#include <Fit/BinData.h>
//#include <Fit/Chi2FCN.h>
//#include <Math/WrappedMultiTF1.h>
//#include <HFitInterface.h>

#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/Algorithms.h"
using namespace RooFit;

struct Result {
	vector<double> Mean_Sigma;
	vector<double> Mean_SquareResidual;
} ;

void PSRanalysis(	TString directory,      // Directory to the data
			TString ConfFile,	// Configuration files
			double Nfilter,		// Filter for the running sum
			double fraction,	// Parameter of Constant Fraction
			double Nsigma,		// Paratemer for t-student test
			double Nthr,		// Parameter of Threshold algorithm
			double thr1,            // first threshold for forward/reversed
			double thr2,            // second threshold for forward/reversed
			TString folder		// Where to save the plots
					){
	
	// Read configuration values
	ReadConfFile Params(ConfFile); // Read the values from the configuration file 
	// Print configuration values
	//Params.Print();
	
	// Compute the baseline
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Baseline
	
	// Dataset List
	vector<string> File = {"scanMvaData3/run1_1.root", "scanMvaData3/run2_1.root"};
	
	// Define some vectors to store the results
	vector<double> MCtruth;             // MC generated values
	vector<double> v1_cfrac, v2_cfrac;  // save constant fraction results
	vector<double> v1_sign, v2_sign;    // save significance results
	
	vector<double> cb_time;             // time of the cb repetitions
	vector<double> da_time;             // time of the da repetitions
	
	
	// REPRODUCE PSR analysis
	ROOT::RDataFrame frame("myTree", {File[0], File[1]}); // Load dataset
	
	for(int i = 0; i < Params.Repetition ; ++i){
	
		TString repetition = TString::Format("repetition == %d", i); // string for Selecting with RdataFrame the repetition
		
		auto timecb = frame.Filter(repetition.Data())
		              .Filter("mwfrequence <= 1000")   // filter on the micro wave frequency to identify the c to b transition
		              .Take<double>("deltaT");  // extract start time of cb transition
		
		auto timeda = frame.Filter(repetition.Data())
		              .Filter("mwfrequence >= 1000")   // filter on the micro wave frequency to identify the d to a transition
		              .Take<double>("deltaT");  // extract start time of da transition
		
		// C TO B ANALYSIS
		auto actualtime_cb = timecb.GetValue(); // get time from Rdataframe node
		cb_time.push_back(actualtime_cb[0]);    // save time in the vector
		
		// Extract from the dataset the starting frequency of the transition
		auto sweepStart_cb = frame.Filter(repetition.Data())
		                     .Filter("mwfrequence <= 1000")
		                     .Take<double>("sweepStart");
		auto actualStart_cb = sweepStart_cb.GetValue(); // get starting frequency from Rdataframe node
		
		// Load the data and save it in a Histogram (raw lineshape)
		 // the histogram start from actualStart_cb, with a step given by Params.FrequencyStep
		auto Spectra1 = frame.Filter(repetition.Data())
		     .Filter("mwfrequence <= 1000")    // filter on the micro wave frequency to identify the c to b transition
		     .Histo1D({"Counts","Frequence", static_cast<int>(Params.SweepStep),actualStart_cb[0], actualStart_cb[0] +  Params.SweepStep* Params.FrequencyStep}, "mwfrequence");

		// ONSET FINDING
		  // Constant fraction algorithm
		double onset1;              // variable where to save the onset
		onset1 = constFrac(Spectra1, fraction, CosmicBackground,Nfilter); 
		v1_cfrac.push_back(onset1); // save the onset in a vector
		
		 // significance algorithm
		onset1 = Significance(Spectra1, Nsigma, CosmicBackground, Nfilter);
		v1_sign.push_back(onset1); // save the onset in a vector
		
		// D TO A ANALYSIS
		auto actualtime_da = timeda.GetValue();  // get time from Rdataframe node
		da_time.push_back(actualtime_da[0]);     // save time in the vector
		
		// Extract from the dataset the starting frequency of the transition	
		auto sweepStart_da = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("sweepStart");
		auto actualStart_da = sweepStart_da.GetValue(); // get starting frequency from Rdataframe node
		
		// Load the data and save it in a Histogram (raw lineshape)
		  // the histogram start from actualStart_da, with a step given by Params.FrequencyStep
		auto Spectra2 = frame.Filter(repetition.Data())
		     .Filter("mwfrequence >= 1000")   // filter on the micro wave frequency to identify the d to a transition
		     .Histo1D({"Counts","Frequence", static_cast<int>( Params.SweepStep), actualStart_da[0], actualStart_da[0] +  Params.SweepStep* Params.FrequencyStep}, "mwfrequence");
		
		// ONSET FINDING
		  // Constant fraction algorithm
		double onset2;	// variable where to save the onset
		onset2 = constFrac(Spectra2, fraction, CosmicBackground,Nfilter);
		v2_cfrac.push_back(onset2);  // save the onset in a vector

		// ONSET FINDING
		  // Significance algorithm
		onset2 = Significance(Spectra2, Nsigma, CosmicBackground, Nfilter);
		v2_sign.push_back(onset2);
		
		// LOAD THE TRUE ONSET VALUE FROM THE DATA
		auto frame3 = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("trueOnset");
		auto frame4 = frame.Filter(repetition.Data()).Filter("mwfrequence <= 1000").Take<double>("trueOnset");
		auto true_da = frame3.GetValue();
		auto true_cb = frame4.GetValue();
		
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back(true_da[0] - true_cb[0]);
		}

	// Combined Fit
   	
   	auto canvas = new TCanvas("canvas", "Constant Fraction", 1000,550);
	auto pad = new TPad("pad", "pad",0,0,1,1);
	
	TMultiGraph *mg = new TMultiGraph();
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad->Draw();
	
	// C to B transition
	auto g = new TGraph(cb_time.size(), cb_time.data(), v1_cfrac.data());
	g->SetMarkerStyle(21); // Medium Dot
	g->SetMarkerColor(9);
	g->SetMarkerSize(1);
	g->SetLineStyle(0);
	// D to A transition
	auto g2 = new TGraph(da_time.size(), da_time.data(), v2_cfrac.data());
	g2->SetMarkerStyle(21); // Medium Dot
	g2->SetMarkerColor(2);
	g2-> SetMarkerSize(1);
	g2->SetLineStyle(0);
	
	mg->SetTitle("Onset vs. Time; time [seconds] ; frequencies ");
	mg->Add(g);
	mg->Add(g2);
	mg->Draw("AP");
	canvas->SaveAs("PSR.pdf");
}

