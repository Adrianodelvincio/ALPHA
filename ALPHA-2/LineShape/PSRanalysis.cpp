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
};


void PSRanalysis(	TString directory,      // Directory to the data
			TString ConfFile,	// Configuration files
			int start,
			int stop,
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
	
	std::vector<std::string> FileList;
	FileList = getFiles(start,stop, directory); // file list to be analyzed
	
	// Compute the baseline
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Baseline
	
	// Dataset List
	//vector<string> File = {"scanMvaData3/run1_1.root", "scanMvaData3/run2_1.root"};
	
	// Define some vectors to store the results
	vector<double> v1_cfrac(Params.Repetition), v2_cfrac(Params.Repetition);  // save constant fraction results
	vector<double> v1_sign(Params.Repetition), v2_sign(Params.Repetition);    // save significance results
	
	vector<double> cb_time(Params.Repetition);             // time of the cb repetitions
	vector<double> da_time(Params.Repetition);             // time of the da repetitions
	
	// Save the results of the fit
	vector<double> vm,vq1,vq2, vdiff;
	vector<double> MCtruth;   // MC generated values
	
	
	// LOOP ON TRIALS
	int trial = 0;
	for(int j = 0; j < FileList.size(); j += 2){
		trial  += 1; std::cout << "Analizzo DataFrame " << trial << "\n" << std::endl;
		
		// REPRODUCE PSR analysis
		ROOT::RDataFrame frame("myTree", {FileList[j], FileList[j+1]}); // Load dataset
		
		// LOOP ON REPETITION
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
			cb_time[i] = (actualtime_cb[0]);    // save time in the vector
			
			// Extract from the dataset the starting frequency of the transition
			auto sweepStart_cb = frame.Filter(repetition.Data())
				             .Filter("mwfrequence <= 1000")
				             .Take<double>("sweepStart");
			auto actualStart_cb = sweepStart_cb.GetValue(); // get starting frequency from Rdataframe node
			
			// Load the data and save it in a Histogram (raw lineshape)
			 // the histogram start from actualStart_cb, with a step given by Params.FrequencyStep
			auto Spectra1 = frame.Filter(repetition.Data())
			     .Filter("mwfrequence <= 1000")    // filter on the micro wave frequency to identify the c to b transition
			     .Histo1D({"Counts","Frequence", static_cast<int>(Params.SweepStep),
			                                     actualStart_cb[0],
			                                     actualStart_cb[0] +  Params.SweepStep* Params.FrequencyStep},
			                                     "mwfrequence");

			// ONSET FINDING
			  // Constant fraction algorithm
			double onset1;              // variable where to save the onset
			onset1 = constFrac(Spectra1, fraction, CosmicBackground,Nfilter); 
			v1_cfrac[i] = (onset1); // save the onset in a vector
			
			 // significance algorithm
			onset1 = Significance(Spectra1, Nsigma, CosmicBackground, Nfilter);
			v1_sign[i] = (onset1); // save the onset in a vector
			
			// D TO A ANALYSIS
			auto actualtime_da = timeda.GetValue();  // get time from Rdataframe node
			da_time[i] = (actualtime_da[0]);     // save time in the vector
			
			// Extract from the dataset the starting frequency of the transition	
			auto sweepStart_da = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("sweepStart");
			auto actualStart_da = sweepStart_da.GetValue(); // get starting frequency from Rdataframe node
			
			// Load the data and save it in a Histogram (raw lineshape)
			  // the histogram start from actualStart_da, with a step given by Params.FrequencyStep
			auto Spectra2 = frame.Filter(repetition.Data())
			     .Filter("mwfrequence >= 1000")   // filter on the micro wave frequency to identify the d to a transition
			     .Histo1D({"Counts","Frequence", static_cast<int>( Params.SweepStep),
			                                     actualStart_da[0],
			                                     actualStart_da[0] +  Params.SweepStep* Params.FrequencyStep},
			                                     "mwfrequence");
			
			// ONSET FINDING
			  // Constant fraction algorithm
			double onset2;	// variable where to save the onset
			onset2 = constFrac(Spectra2, fraction, CosmicBackground,Nfilter);
			v2_cfrac[i] = (onset2);  // save the onset in a vector

			// ONSET FINDING
			  // Significance algorithm
			onset2 = Significance(Spectra2, Nsigma, CosmicBackground, Nfilter);
			v2_sign[i] = (onset2);
		} // LOOP ON REPETITION
		
		/////////////////////////////////////////
		// Combined Fit
	   	double sigmay = 3.089;
	   	//fit to find the slope
	   	double slope = 0;
	   	double S1xy = 0;
	   	double S0x = 0;
	   	double S0xy = 0;
	   	double S1x = 0;
	   	double S2x = 0;
	   	
	   	double meancbtime = 0;
	   	double meandatime = 0;
	   	for(int i = 0; i < Params.Repetition; i++){
	   		meancbtime += (cb_time[i])/Params.Repetition;
	   		meandatime += (da_time[i])/Params.Repetition;
	   	}
	   	double difftime = meancbtime - meandatime;
	   	for(int i = 0; i < Params.Repetition; i++){
	   		// the quantities are computer with an allinement, not formal !!!
	   		S1xy += (1/pow(sigmay,2)) * (cb_time[i]*v1_cfrac[i] + (da_time[i] + difftime)*v2_cfrac[i]);
	   		S0x  += 1/pow(sigmay,2) + 1/pow(sigmay,2);
	   		S0xy += 1/pow(sigmay,2) * (v1_cfrac[i] + v2_cfrac[i]);
	   		S1x  += 1/pow(sigmay,2) * (cb_time[i] + (da_time[i] + difftime));
	   		S2x  += 1/pow(sigmay,2) * (cb_time[i]*cb_time[i] + (da_time[i] + difftime)*(da_time[i] + difftime));
	   	}

	   	double D = S0x*S2x - S1x*S1x;
	   	slope = (S1xy*S0x - S0xy*S1x)/D;
	   	double dslope = S0x/D;
	   	std::cout << "m is: " << slope << " +/- " << sqrt(dslope) << std::endl;
		
		// fit for the onset
	   	TF1 *f1 = new TF1("f1", "[0] + [1]*x",0,1e6);
		f1->FixParameter(1,slope);
		TF1 *f2 = new TF1("f2", "[0] + [1]*x",0,1e6);
		f2->FixParameter(1,slope);
		
		f1->SetParName(0, "q1");
		f2->SetParName(0, "q2");
		f1->SetParName(1, "m");
		f2->SetParName(1, "m");
		
		vector<double> ex(Params.Repetition,0);
	   	vector<double> ey(Params.Repetition,sigmay);
	   	
	   	auto g = new TGraphErrors(cb_time.size(), cb_time.data(), v1_cfrac.data(), ex.data(), ey.data());
	   	g->Fit(f1,"","same");
		TF1 *fit1 = g->GetFunction("f1");
		double q1 = fit1->GetParameter(0);
		double eq1 = fit1->GetParError(0);
	   	auto g2 = new TGraphErrors(da_time.size(), da_time.data(), v2_cfrac.data(), ex.data(), ey.data());
	   	g2->Fit(f2,"","same");
		TF1 *fit2 = g2->GetFunction("f2");
		//double chi2_2 = fit2->GetChisquare();
		double q2 = fit2->GetParameter(0);
		double eq2 = fit2->GetParError(0);
		
		std::setprecision(8);
		std::cout << "intercept cb: " << q1 << " +/- " << eq1 << std::endl;
		std::cout << "intercept da: " << q2 << " +/- " << eq2 << std::endl;
		std::cout << "hyperfine measured : " << std::setprecision(8) << q2 - q1 << " +/- " << eq1 + eq2 << std::endl; 
		/////////////////////////////////////////
		
		// save the data
		
		// LOAD THE TRUE ONSET VALUE FROM THE DATA
		TString repetition = TString::Format("repetition == 0"); // string for Selecting with RdataFrame the repetition
		auto frame3 = frame.Filter(repetition.Data()).Filter("mwfrequence >= 1000").Take<double>("trueOnset");
		auto frame4 = frame.Filter(repetition.Data()).Filter("mwfrequence <= 1000").Take<double>("trueOnset");
		auto true_da = frame3.GetValue();
		auto true_cb = frame4.GetValue();
		
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back(true_da[0] - true_cb[0]);
		vm.push_back(slope);
		vq1.push_back(q1);
		vq2.push_back(q2);
		vdiff.push_back((q2 - q1) - (true_da[0] - true_cb[0]));
	} // LOOP ON TRIALS
	
	// Montecarlo
	std::vector<double> w(vm.size(),1); // weights vector
	
	auto h1 = new TH1D("h1","q1", 81, 0, 100);
	auto h2 = new TH1D("h2","q2", 81, 1420000,1420000 + 100);
	auto h3 = new TH1D("h3","slope",1000, -10 - 5e-4 , - 10 + 5e-4);
   	auto h4 = new TH1D("h3","hfs_{measured} - hfs_{generated}", 41, -20 , 20);
   	
   	h1->FillN(vq1.size(),vq1.data(), w.data());
   	h2->FillN(vq2.size(),vq2.data(), w.data());
   	h3->FillN(vm.size(),vm.data(), w.data());
	h4->FillN(vm.size(),vdiff.data(), w.data());
	
	auto Montecarlo = new TCanvas("montecarlo", "Montecarlo Fit", 1000,550);
	auto pad2 = new TPad("pad2", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad2->Divide(2,2,0.001,0.001);
	pad2->Draw();
	pad2->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
	pad2->cd(2);
	h2->GetXaxis()->SetTitle("frequency [kHz]");
	h2->SetLineColor(2);
	h2->SetLineWidth(2);
	h2->Draw();
	pad2->cd(3);
	h3->GetXaxis()->SetTitle("frequency/second [kHz/second]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	pad2->cd(4);
	h4->GetXaxis()->SetTitle("frequency [kHz]");
	h4->SetLineWidth(2);
	h4->SetLineColor(38);
	h4->Draw();
	
	
	
	// Combined Fit
   	double sigmay = 3.089;
   	vector<double> ex(Params.Repetition,0);
	vector<double> ey(Params.Repetition,sigmay);
   	//fit to find the slope
   	double slope = 0;
   	double S1xy = 0;
   	double S0x = 0;
   	double S0xy = 0;
   	double S1x = 0;
   	double S2x = 0;

   	double meancbtime = 0;
   	double meandatime = 0;
   	for(int i = 0; i < Params.Repetition; i++){
   		meancbtime += (cb_time[i])/Params.Repetition;
   		meandatime += (da_time[i])/Params.Repetition;
   	}
   	double difftime = meancbtime - meandatime;
   	for(int i = 0; i < Params.Repetition; i++){
   		// the quantities are computer with an allinement, not formal !!!
   		S1xy += (1/pow(sigmay,2)) * (cb_time[i]*v1_cfrac[i] + (da_time[i] + difftime)*v2_cfrac[i]);
   		S0x  += 1/pow(sigmay,2) + 1/pow(sigmay,2);
   		S0xy += 1/pow(sigmay,2) * (v1_cfrac[i] + v2_cfrac[i]);
   		S1x  += 1/pow(sigmay,2) * (cb_time[i] + (da_time[i] + difftime));
   		S2x  += 1/pow(sigmay,2) * (cb_time[i]*cb_time[i] + (da_time[i] + difftime)*(da_time[i] + difftime));
   	}

   	double D = S0x*S2x - S1x*S1x;
   	slope = (S1xy*S0x - S0xy*S1x)/D;
   	double dslope = S0x/D;
   	std::cout << "m is: " << slope << " +/- " << sqrt(dslope) << std::endl;

	// fit for the onset
   	TF1 *f1 = new TF1("f1", "[0] + [1]*x",0,1e6);
	f1->FixParameter(1,slope);
	TF1 *f2 = new TF1("f2", "[0] + [1]*x",0,1e6);
	f2->FixParameter(1,slope);
	
	f1->SetParName(0, "q1");
	f2->SetParName(0, "q2");
	f1->SetParName(1, "m");
	f2->SetParName(1, "m");

   	auto canvas = new TCanvas("canvas", "Constant Fraction", 1000,550);
	auto pad = new TPad("pad", "pad",0,0,1,1);
	gStyle->SetOptFit(1);
	TMultiGraph *mg = new TMultiGraph();
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad->Draw();
	
	// C to B transition
	auto g = new TGraphErrors(cb_time.size(), cb_time.data(), v1_cfrac.data(), ex.data(), ey.data());
	g->SetMarkerStyle(21); // Medium Dot
	g->SetMarkerColor(9);
	g->SetMarkerSize(1);
	g->SetLineStyle(0);
	g->Fit(f1,"+","same");
	TF1 *fit1 = g->GetFunction("f1");
	//double chi2 = fit1->GetChisquare();
	double q1 = fit1->GetParameter(0);
	double eq1 = fit1->GetParError(0);

	// D to A transition
	auto g2 = new TGraphErrors(da_time.size(), da_time.data(), v2_cfrac.data(), ex.data(), ey.data());
	g2->SetMarkerStyle(21); // Medium Dot
	g2->SetMarkerColor(2);
	g2->SetMarkerSize(1);
	g2->SetLineStyle(0);
	g2->Fit(f2,"+","same");
	TF1 *fit2 = g2->GetFunction("f2");
	//double chi2_2 = fit2->GetChisquare();
	double q2 = fit2->GetParameter(0);
	double eq2 = fit2->GetParError(0);
	
	std::setprecision(8);
	std::cout << "intercept cb: " << q1 << " +/- " << eq1 << std::endl;
	std::cout << "intercept da: " << q2 << " +/- " << eq2 << std::endl;
	std::cout << "hyperfine measured : " << std::setprecision(8) << q2 - q1 << " +/- " << eq1 + eq2 << std::endl; 
	
	mg->SetTitle("Onset vs. Time; time [seconds] ; frequencies ");
	mg->Add(g);
	mg->Add(g2);
	mg->Draw("AP");
	canvas->SaveAs("PSR.pdf");	
}

