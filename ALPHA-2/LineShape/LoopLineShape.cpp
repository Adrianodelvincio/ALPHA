#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TSpline.h"
#include <TMath.h>
#include <TRandom3.h>
#include <math.h> 
#include <cmath>
#include "../Headers/toyLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/LineShape.h"

using namespace RooFit;

void LoopLineShape(	int Nloop,
			TString folder,
			TString mvaScan,
			double CosmicRate,
			double Efficiency,
			TString ConfFile){

// ReadConfigurationFiles;
ReadConfFile Params(ConfFile);

// Select the efficiency directly form args of the macro
if(mvaScan == TString::Format("true")){
	Params.Efficiency = Efficiency;
	Params.CosmicRate = CosmicRate;
}
//Params.Print();


// Compute expected number of cosmic annihilation and background
double Ntot = Params.Nstack * Params.NHbar * Params.Efficiency;		// Number of Total Events
double Nc = Ntot*Params.C; 				// Expected event for lineshape1
double Nd = Ntot*(1 - Params.C); 			// Expected event for lineshape2
double Ncosmic = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
//	Lineshape Parameters
double cb_start = Params.cb_start;
double cb_end = Params.cb_end;
double cb_peak = Params.cb_peak;
double da_start = Params.da_start;
double da_end = Params.da_end;
double da_peak = Params.da_peak;
//	Cruijff Parameters
double k0_cb = Params.k0_cb;
double k1_cb = Params.k1_cb;
double sigma0_cb = Params.sigma0_cb;
double sigma1_cb = Params.sigma1_cb;
double Norm_cb = Params.Norm_cb;
double k0_da = Params.k0_da;
double k1_da = Params.k1_da;
double sigma0_da = Params.sigma0_da;
double sigma1_da = Params.sigma1_da;
double Norm_da = Params.Norm_da;


int Ntrial = Nloop;				// Ntrial
double pGas_d = 1 - Params.WallComponent_ad; 			// Percentage of annihilation on residual gas
double pGas_c = 1 - Params.WallComponent_cb; 			// Percentage of annihilation on residual gas

int Nbin1 = Params.TotalStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP
int Nbin2 = Params.TotalStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

RooRealVar x("x","r [cm]",0.,4.);
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85, 0.85, 0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix); //Pdf Annihilation on walls
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay)); //PDF Residual Gas
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x)); //PDF Cosmic

RooRealVar Nwall("Nwall","Nwall",100, -3000, +1000000);
RooRealVar Ngas ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk 	("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model","model", RooArgList{gauss_Mix}, RooArgList{Nwall});
RooAddPdf genGas("model1","model1", RooArgList{Rayleigh}, RooArgList{Ngas});
RooAddPdf genCosmic("model2", "model2", RooArgList{linearFit}, RooArgList{Nbk});

/////////////////
// EXTERNAL TOY LOOP
TRandom3 *r = new TRandom3();
	for(int l = 0; l < Ntrial; l++){				//LOOP ON TRIALS
		
		// INITIALIZATION OF VARIABLES
		RooDataSet dataPdf1("data", "data", RooArgSet(x));	// Dataset to store the events
		RooDataSet dataPdf2("dati", "dati", RooArgSet(x));	// Dataset to store the events
		vector<double>	f1, f2;					// Frequencies pdf1 , Frequencies pdf2
		vector<double>	v1Tot, v2Tot;				// Total Counts
		vector<int>	v1Type,v2Type;				// Type of Event
		vector<int>	RunNumber1, RunNumber2;			// Run Number (from 0 to Repetition)
		vector<double>  generated_cb;				// shift of c to b onset
		vector<double>	generated_da;				// shift of d to a onset
		vector<double>  start_sweep_cb;
		vector<double>  start_sweep_da;
		vector<double>  time_cb;
		vector<double>  time_da;
		
		// HYPERFINE SPLITTING GENERATION
		double delta = r->Uniform(Params.delta_left,Params.delta_right); 		// generate a shift in the hyperfine splitting
		double LineShift_cb = -(delta/2);	// shift of the c to b lineshape
		double LineShift_da = +(delta/2);	// shift of the d to a lineshape
		//
		
		// OPERATOR SET THE WINDOW BEFORE THE FIRST RUN, DEFINITION OF THE TRUE Bdrift
		double OperatorError =  r->Uniform(Params.FrequencyMinBeforeOnset,Params.FrequencyMaxBeforeOnset);
		double TrueBdrift = Params.Bdrift + r->Gaus(0,Params.sigmaDrift);
		//
		
		// Define time 
		double TimeInBetween = 0; // In seconds
		//
		
		// LOOP ON REPETITION
		for(int run = 0; run < Params.Repetition; run++){
			//////////////
			// C TO B TRANSITION
			
			// save time
			time_cb.push_back(TimeInBetween);
			
			// generate and save true onset
			double onset_cb = cb_start + LineShift_cb - TrueBdrift*(TimeInBetween); 	// onset c to b + shift + drift
			double cbMaximum  = cb_peak  + LineShift_cb - TrueBdrift*(TimeInBetween);	// peak c to b + shift + drift	
			//generated_cb.push_back(cb_start + LineShift_cb - TrueBdrift*(TimeInBetween)); 	
			generated_cb.push_back(cb_start + LineShift_cb); 
			
			// Start of frequency sweep c-b			
			double startSweepCB = Params.cb_start + OperatorError - Params.Bdrift*(TimeInBetween);
			start_sweep_cb.push_back(startSweepCB + Params.FrequencyStep/2);
			
			// Create the series of micro wave measurements for c to b
			TH1F *genLineShape1 = new TH1F("hist1", "lineshape c to b", Nbin1, startSweepCB, startSweepCB + Nbin1*(Params.FrequencyStep));
			SetContent(genLineShape1,Nbin1,Cruijff, onset_cb, cbMaximum, sigma0_cb, sigma1_cb, k0_cb, k1_cb, Norm_cb);

			// update time
			TimeInBetween += (Params.SweepStep + Params.ClearingStep)*Params.TimeStep; 
			
			//////////////
			// D TO A TRANSITION
			// save time
			time_da.push_back(TimeInBetween);
			
			// generate and save true onset
			double onset_da = da_start + LineShift_da - TrueBdrift*(TimeInBetween); 	// onset d to a + shift + drift
			double daMaximum  = da_peak  + LineShift_da - TrueBdrift*(TimeInBetween); 	// peak d to a + shift + drift
			//generated_da.push_back(da_start + LineShift_da - TrueBdrift*(TimeInBetween));
			generated_da.push_back(da_start + LineShift_da);
			
			// Start of frequency sweep d-a
			double startSweepDA = Params.da_start + OperatorError - Params.Bdrift*(TimeInBetween);
			start_sweep_da.push_back(startSweepDA + Params.FrequencyStep/2);

			// Create the series of micro wave measurements for d to a
			TH1F *genLineShape2 = new TH1F("hist2", "lineshape d to a", Nbin2, startSweepDA, startSweepDA + Nbin2*(Params.FrequencyStep)); // create the transition
			SetContent(genLineShape2,Nbin2,Cruijff, onset_da, daMaximum, sigma0_da, sigma1_da, k0_da, k1_da, Norm_da);

			// GENERATE THE DATA, LOOP ON BINS
			// time relative to the sweep start for each bin
			double tbc = 0;
			double tda = 0;
			
			for(int bin = 1; bin <= Params.SweepStep; bin++){
				
				// Bdrift correction for each bin
				//onset_cb  += - TrueBdrift*(Params.TimeStep);  // onset c to b + shift + drift
				//cbMaximum += - TrueBdrift*(Params.TimeStep);	// peak c to b + shift + drift	
				//onset_da  += - TrueBdrift*(Params.TimeStep); 	// onset d to a + shift + drift
				//daMaximum += - TrueBdrift*(Params.TimeStep); 	// peak d to a + shift + drift
				
				// update the lineshape
				SetContent(genLineShape1,Nbin1,Cruijff, onset_cb, cbMaximum, sigma0_cb, sigma1_cb, k0_cb, k1_cb, Norm_cb);
				SetContent(genLineShape2,Nbin2,Cruijff, onset_da, daMaximum, sigma0_da, sigma1_da, k0_da, k1_da, Norm_da);
				
				/////////
				// LINESHAPE C TO B
				
				// ANNIHILATION ON WALLS
				double prob = ComputeProb(genLineShape1,bin);// Probability of the bin
				SetCoefficients((Params.WallComponent_cb*Nc)*prob,(pGas_c*Nc)/Params.SweepStep, Ncosmic, &Nwall,&Ngas,&Nbk);
				int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
				double frequence = genLineShape1->GetBinCenter(bin);
				
				if(prob > 0){
				RooDataSet *dataLoopWall = genMix.generate(x,Extended()); // GENERATE THE DATA
				SetVectors(dataPdf1,dataLoopWall,
							v1Type, 0,
							mixCount,
							f1,frequence,
							RunNumber1,run); 
				delete dataLoopWall;
				}else { mixCount = 0;}
				// ANNIHILATION DUE TO RESIDUAL GAS
				if(Ngas.getVal() != 0){
				RooDataSet *dataLoopGas = genGas.generate(x, Extended());
				SetVectors(dataPdf1,dataLoopGas,
							v1Type, 1,
							gasCount,
							f1,frequence,
							RunNumber1,run);
				delete dataLoopGas;
				}else {gasCount = 0;}
				// COSMIC BACKGROUND
				if(Ncosmic > 0){
				RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
				SetVectors(dataPdf1,dataCosmic,
							v1Type, 2,
							CosmicCount,
							f1,frequence,
							RunNumber1,run);
				delete dataCosmic;
				}else{ CosmicCount = 0;}
				// STORE SOME USEFUL QUANTITIES
				v1Tot.push_back(mixCount + gasCount + CosmicCount);
				//std::cout << "line1 " << "Bin: " << bin << " frequence: " << frequence << " " << mixCount << " " << CosmicCount << std::endl;
				
				/////////////////
				// LINESHAPE D TO A
				mixCount = 0; gasCount = 0; CosmicCount = 0;
				// ANNIHILATION ON WALLS
				prob = ComputeProb(genLineShape2,bin);	// Probability of the bin
				frequence = genLineShape2->GetBinCenter(bin);
				SetCoefficients((Params.WallComponent_ad*Nd)*prob,(pGas_d*Nd)/Params.SweepStep,Ncosmic, &Nwall,&Ngas,&Nbk);
				
				if(prob > 0){
				RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
				SetVectors(dataPdf2,dataLoopWall2,
							v2Type,0,
							mixCount,
							f2,frequence,
							RunNumber2,run);
				delete dataLoopWall2;
				} else{  mixCount = 0;}
				// ANNIHILATION DUE TO RESIDUAL GAS
				if(Ngas.getVal() != 0){
				RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
				SetVectors(dataPdf2,dataLoopGas2,
							v2Type,1,
							gasCount,
							f2,frequence,
							RunNumber2,run);
				delete dataLoopGas2;
				} else{ gasCount = 0;}
				// COSMIC BACKGROUND
				if(Ncosmic > 0){
				RooDataSet *dataCosmic2 = genCosmic.generate(x, Extended());			
				SetVectors(dataPdf2,dataCosmic2,
							v2Type,2,
							CosmicCount,
							f2,frequence,
							RunNumber2, run);
				delete dataCosmic2;
				}else{ CosmicCount = 0;}
				// STORE USEFUL QUANTITIES
				v2Tot.push_back(mixCount + gasCount + CosmicCount);
				//std::cout << "line2 " << "Bin: " << bin << " frequence: " << frequence - 1420000 << " " << mixCount << " " << CosmicCount << std::endl;
				
				// update time
				tda += Params.TimeStep;
				tbc += Params.TimeStep;
			} // Loop on Bin
			
			//update time for new repetition
			TimeInBetween += (1.30 + r->Uniform(0,0.45))*3600;
			
			delete genLineShape1;
			delete genLineShape2;
		} // Loop on Repetition
	
		int Tot1 = std::accumulate(v1Tot.begin(), v1Tot.end(), 0);
		int Tot2 = std::accumulate(v2Tot.begin(), v2Tot.end(), 0);
		// SAVE SIMULATED EVENTS 
		ROOT::RDataFrame d1(Tot1 -1);	// PDF1
		ROOT::RDataFrame d2(Tot2 -1); 	// PDF2
		
		TString nameFile1 = TString::Format("run1_%d.root", l);
		TString nameFile2 = TString::Format("run2_%d.root", l);
		
		int j(0);
		auto FilledFrame1 = FillDataFrame(d1,
						dataPdf1,
						f1,
						v1Type,
						v1Tot,
						j,
						Params.TimeStep,
						Params.FrequencyStep,
						start_sweep_cb,
						time_cb,
						RunNumber1,
						generated_cb); // Fill RDataFrame
		std::cout << "Save " <<  folder + nameFile1 << std::endl;
		FilledFrame1.Snapshot("myTree", folder + nameFile1);

		j = 0; // FROM HERE FILL DA DATA
		auto FilledFrame2 = FillDataFrame(d2,
						dataPdf2,
						f2,
						v2Type,
						v2Tot,
						j,
						Params.TimeStep,
						Params.FrequencyStep,
						start_sweep_da,
						time_da,
						RunNumber2,
						generated_da); // Fill RDataFrame
		std::cout << "Save " << folder + nameFile2 << std::endl;
		FilledFrame2.Snapshot("myTree", folder + nameFile2);
		//std::cout << "Efficiency: " << Params.Efficiency << " Background: " << Params.CosmicRate << std::endl;
		//std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
		//std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl;
		dataPdf1.Delete(); dataPdf2.Delete();
	} // Loop Trials
	
	// Show the lineshape
	/*
	SetContent(genLineShape1,Nbin1,Cruijff, cb_start, cb_peak, sigma0_cb, sigma1_cb, k0_cb, k1_cb, Norm_cb);
	SetContent(genLineShape2,Nbin2,Cruijff, da_start, da_peak, sigma0_da, sigma1_da, k0_da, k1_da, Norm_da);
	TH1F *genLineShape3 = new TH1F("hist3", "lineshape", 500,0,400);
	TH1F *genLineShape4 = new TH1F("hist4", "lineshape", 500,0,400);
	SetContent(genLineShape3,500,Cruijff, 0, 220, sigma0_cb, sigma1_cb, k0_cb, k1_cb, Norm_cb);
	SetContent(genLineShape4,500,Cruijff, 10, 220+10 , sigma0_da, sigma1_da, k0_da, k1_da, Norm_da);
	auto b = new TCanvas("b1", "Spectral lines");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	gStyle->SetOptStat(0);
	genLineShape1->GetXaxis()->SetTitle("frequency [kHz]");
	genLineShape1->SetMarkerStyle(21);
	genLineShape1->SetMarkerColor(1);
	genLineShape1->SetMarkerSize(0.5);
	genLineShape1->SetLineColor(1);
	genLineShape1->Draw();
	//genLineShape3->SetMarkerStyle(21);
	//genLineShape3->SetMarkerColor(kViolet);
	//genLineShape3->SetLineColor(kViolet);
	//genLineShape3->SetMarkerSize(0.5);
	//genLineShape3->Draw("samehist");
	pad->cd(2);
	genLineShape2->GetXaxis()->SetTitle("frequency [kHz]");
	genLineShape2->SetMarkerStyle(21);
	genLineShape2->SetMarkerColor(1);
	genLineShape2->SetMarkerSize(0.5);
	genLineShape2->SetLineColor(1);
	genLineShape2->Draw();
	//genLineShape4->SetMarkerStyle(21);
	//genLineShape4->SetMarkerColor(kViolet);
	//genLineShape4->SetMarkerSize(0.5);
	//genLineShape4->SetLineColor(kViolet);
	//genLineShape4->Draw("samehist");
	auto c = new TCanvas("c1","simone is runnuning me out of patience");
	genLineShape3->Draw("histl");
	genLineShape4->Draw("samehistl");
	*/
}
