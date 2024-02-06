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

void PSRanalysis(	TString directory,
					TString ConfFile,	// Configuration files
					double Nfilter,		// Filter for the running sum
					double fraction,	// Parameter of Constant Fraction
					double Nsigma,		// Paratemer for t-student test
					double Nthr,		// Parameter of Threshold algorithm
					double thr1,
					double thr2,
					TString folder		// Where to save the plots
					){
	
	
	ReadConfFile Params(ConfFile); // Read the values from the configuration file 
	//Params.Print();
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
	double FrequencyStep = Params.FrequencyStep;
	double startPdf1 = Params.x_cb_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.x_da_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep d-a
	double SweepStep = Params.SweepStep; // number of bin for each lineshape
	
	vector<string> File = {"scanMvaData3/run1_1.root", "scanMvaData3/run2_1.root"};
	// Define some vectors to store the results
	vector<double> MCtruth;
	vector<double> v1_cfrac, v2_cfrac;
	vector<double> onsets;
	vector<double> v1_sign, v2_sign;

	vector<double> time;
	vector<double> cb_time;
	vector<double> da_time;
	double threshold = Nthr*CosmicBackground; // threshold considering the cosmic background
	
	// Show the first two files
	ROOT::RDataFrame rdf("myTree", {File[0], File[1]});
	auto hist_ctob = rdf.Filter("repetition == 0")
			.Filter("mwfrequence <= 1000")
			.Histo1D({"Counts"," c to b",static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
	
	auto histF4 = 	rdf.Filter("repetition == 0")
			.Filter("mwfrequence >= 1000")
			.Histo1D({"Counts"," d to a",static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");
	
	auto b = new TCanvas("b1", "Spectral lines");
	auto pad1 = new TPad("pad1", "pad",0,0,1,1);
	pad1->Divide(2,1,0.001,0.001); pad1->Draw();
	pad1->cd(1);
	hist_ctob->SetTitle("c to b without cosmic events");
	hist_ctob->SetMarkerStyle(21);
	hist_ctob->SetMarkerColor(2);
	hist_ctob->SetLineColor(4);
	hist_ctob->DrawClone();
	hist_ctob->DrawClone("SAME P");
	pad1->cd(2);
	histF4->SetTitle("d to a with cosmic events");
	histF4->SetMarkerStyle(21);
	histF4->SetMarkerColor(2);
	histF4->SetLineColor(4);
	histF4->DrawClone();
	histF4->DrawClone("SAME P");
	//
	
	
	// REPRODUCE PSR analysis
	ROOT::RDataFrame frame("myTree", {File[0], File[1]}); // Load dataset
	
	for(int i = 0; i < Params.Repetition ; ++i){
		
		TString repetition = TString::Format("repetition == %d", i); // repetition
		auto timecb = frame.Filter(repetition.Data()).Filter("mwfrequence <= 1000").Take<double>("deltaT");
		auto timeda = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("deltaT");
		// C TO B ANALYSIS
		
		auto actualtime_cb = timecb.GetValue();
		
		time.push_back(actualtime_cb[0]);
		cb_time.push_back(actualtime_cb[0]);
		
		// starting frequency of the sweep 
		auto sweepStart_cb = frame.Filter(repetition.Data()).Filter("mwfrequence <= 1000").Take<double>("sweepStart");
		auto actualStart_cb = sweepStart_cb.GetValue();
		
		// load the sweep in a histogram
		auto Spectra1 = frame.Filter(repetition.Data())
		 	.Filter("frequence <= 1000")
		 	.Histo1D({"Counts","Frequence", static_cast<int>(SweepStep),actualStart_cb[0], actualStart_cb[0] + SweepStep*FrequencyStep}, "frequence");
		
		// measure and save the onset
		double onset1;				// Reconstructed onset
		onset1 = constFrac(Spectra1, fraction, CosmicBackground,Nfilter);
		v1_cfrac.push_back(onset1);
		onsets.push_back(onset1);
		std::cout << actualtime_cb[0] << " onset cb: " << onset1  << std::endl;
		onset1 = Significance(Spectra1, Nsigma, CosmicBackground, Nfilter);
		v1_sign.push_back(onset1);
		
		
		// D TO A ANALYSIS
		auto actualtime_da = timeda.GetValue();
		time.push_back(actualtime_da[0]);
		da_time.push_back(actualtime_da[0]);
		
		// starting frequency of the sweep 		
		auto sweepStart_da = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("sweepStart");
		auto actualStart_da = sweepStart_da.GetValue();
		
		// load the sweep in a histogram		
		auto Spectra2 = frame.Filter(repetition.Data())
		 .Filter("mwfrequence >= 1000")
		 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep), actualStart_da[0], actualStart_da[0] + SweepStep*FrequencyStep}, "frequence");
		
		// measure and save the onset
		double onset2;	// Reconstructed onset
		onset2 = constFrac(Spectra2, fraction, CosmicBackground,Nfilter);
		v2_cfrac.push_back(onset2);
		onsets.push_back(onset2);
		std::cout << actualtime_da[0] << " onset da: " << onset2  << std::endl;
		onset2 = Significance(Spectra2, Nsigma, CosmicBackground, Nfilter);
		v2_sign.push_back(onset2);
		/**/
		// Load the shifts of the lineshapes
		auto frame3 = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("trueOnset");
		auto frame4 = frame.Filter(repetition.Data()).Filter("mwfrequence <= 1000").Take<double>("trueOnset");
		auto true_da = frame3.GetValue();
		auto true_cb = frame4.GetValue();
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back((true_da[0] - true_cb[0]));
		}
   	
   	auto canvas = new TCanvas("canvas", "Constant Fraction", 1000,550);
	auto pad = new TPad("pad", "pad",0,0,1,1);
	
	TMultiGraph *mg = new TMultiGraph();
	
	//canvas->SetLogy();
	//gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	//pad2->Divide(2,2,0.001,0.001);
	pad->Draw();
	auto g = new TGraph(cb_time.size(), cb_time.data(), v1_cfrac.data());
	g->SetMarkerStyle(21); // Medium Dot
	g->SetMarkerColor(9);
	g->SetMarkerSize(1);
	g->SetLineStyle(0);
	
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

