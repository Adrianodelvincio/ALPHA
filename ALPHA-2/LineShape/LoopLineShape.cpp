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
//Params.Print();
if(mvaScan == TString::Format("true")){
	Params.Efficiency = Efficiency;
	Params.CosmicRate = CosmicRate;
}

//	Parameters of the Simulation 
double Ntot = Params.Nstack * Params.NHbar * Params.Efficiency;		// Number of Total Events
double FrequencyStep = Params.FrequencyStep;		// Kilo hertz
int SweepStep = Params.SweepStep;					// Number of FrequencyStep
int BeforeOnset = Params.BinBeforeOnset;
int Repetition = Params.Repetition;					// Repetition of the single run
double pWall_c = Params.WallComponent_cb;			// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Params.WallComponent_ad;			// Weight annihilation on walls for pdf2 (transition d -> a)
double Ncosmic = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
//	Lineshape Parameters
double x_cb_start = Params.x_cb_start;
double x_cb_end = Params.x_cb_end;
double x_cb_peak = Params.x_cb_peak;
double x_da_start = Params.x_da_start;
double x_da_end = Params.x_da_end;
double x_da_peak = Params.x_da_peak;
double range = Params.delay;						// Set range delay
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


int Ntrial = Nloop;						// Ntrial
double Nc = Ntot*Params.C; 				// Expected event for lineshape1
double Nd = Ntot*(1 - Params.C); 		// Expected event for lineshape2
double pGas_d = 1 - pWall_d; 			// Percentage of annihilation on residual gas
double pGas_c = 1 - pWall_c; 			// Percentage of annihilation on residual gas

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


///////////////////
// DISCRETIZED CONTAINER OF THE LINESHAPE
double ScanStart1 = Params.x_cb_start - (FrequencyStep)*(BeforeOnset + 0.5);	// Start of frequency sweep c-b
double ScanStart2 = Params.x_da_start - (FrequencyStep)*(BeforeOnset + 0.5);	// Start of frequency sweep d-a
TH1F *genLineShape1 = new TH1F("hist1", "lineshape c to b", Nbin1, ScanStart1, ScanStart1 + Nbin1*(FrequencyStep));
TH1F *genLineShape2 = new TH1F("hist2", "lineshape d to a", Nbin2, ScanStart2, ScanStart2 + Nbin2*(FrequencyStep));
//////////////////

/////////////////
// EXTERNAL TOY LOOP
TRandom3 *r = new TRandom3();

	for(int l = 0; l < Ntrial; l++){				//LOOP ON TRIALS
		RooDataSet dataPdf1("data", "data", RooArgSet(x));	// Dataset to store the events
		RooDataSet dataPdf2("dati", "dati", RooArgSet(x));	// Dataset to store the events
		vector<double>	f1, f2;					// Frequencies pdf1 , Frequencies pdf2
		vector<double>	v1Tot, v2Tot;			// Total Counts
		vector<int>	v1Type,v2Type;				// Type of Event
		vector<int>	RunNumber1, RunNumber2;		// Run Number (from 0 to Repetition)
		vector<double>  cb_shift;				// shift of c to b onset
		vector<double>	da_shift;				// shift of d to a onset
		
		//double shift = r->Uniform(-3*Params.FrequencyStep,+3*Params.FrequencyStep); // generate a shift in the hyperfine splitting
		double shift = r->Uniform(-6*Params.FrequencyStep,+2*Params.FrequencyStep); // generate a shift in the hyperfine splitting
		double LineShift_cb = -(shift/2);	// shift of the c to b lineshape
		double LineShift_da = +(shift/2);	// shift of the d to a lineshape
		
		//std::cout << "LineShift c to b " << LineShift_cb << std::endl;
		//std::cout << "LineShift d to a " << LineShift_da << std::endl;
		/////////////
		// LOOP ON REPETITION
		for(int run = 0; run < Repetition; run++){
			//double smearing = r->Uniform(-range, range); 	// set smearing of the onset inside frequency step
			cb_shift.push_back(LineShift_cb);				// Save shift
			da_shift.push_back(LineShift_da);				// Save shift
			double start1 = x_cb_start + LineShift_cb;
			double peak1  = x_cb_peak  + LineShift_cb;
			double start2 = x_da_start + LineShift_da;
			double peak2  = x_da_peak  + LineShift_da;
			
			///////////////
			// DISCTRETIZATION OF THE LINESHAPE
			SetContent(genLineShape1,Nbin1,Cruijff, start1, peak1, sigma0_cb, sigma1_cb, k0_cb, k1_cb, Norm_cb);
			SetContent(genLineShape2,Nbin2,Cruijff, start2, peak2, sigma0_da, sigma1_da, k0_da, k1_da, Norm_da);
			
			for(int bin = 1; bin <= SweepStep; bin++){ 	//LOOP ON BINS
				// PDF1
				double prob = ComputeProb(genLineShape1,bin);// Probability of the bin
				SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/SweepStep, Ncosmic, &Nwall,&Ngas,&Nbk);
				int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
				double frequence = genLineShape1->GetBinCenter(bin);
				
				if(prob > 0){
				RooDataSet *dataLoopWall = genMix.generate(x,Extended()); // GENERATE THE DATA
				SetVectors(dataPdf1,dataLoopWall,
							v1Type, 0,
							mixCount,
							f1,frequence,
							RunNumber1,run); // FILL DATASET
				delete dataLoopWall;
				}else { mixCount = 0;}
				if(Ngas.getVal() != 0){
				RooDataSet *dataLoopGas = genGas.generate(x, Extended());
				SetVectors(dataPdf1,dataLoopGas,
							v1Type, 1,
							gasCount,
							f1,frequence,
							RunNumber1,run);
				delete dataLoopGas;
				}else {gasCount = 0;}
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
				// PDF 2
				mixCount = 0; gasCount = 0; CosmicCount = 0;
				prob = ComputeProb(genLineShape2,bin);	// Probability of the bin
				frequence = genLineShape2->GetBinCenter(bin);
				SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/SweepStep,Ncosmic, &Nwall,&Ngas,&Nbk);
				
				if(prob > 0){
				RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
				SetVectors(dataPdf2,dataLoopWall2,
							v2Type,0,
							mixCount,
							f2,frequence,
							RunNumber2,run);
				delete dataLoopWall2;
				} else{  mixCount = 0;}
				if(Ngas.getVal() != 0){
				RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
				SetVectors(dataPdf2,dataLoopGas2,
							v2Type,1,
							gasCount,
							f2,frequence,
							RunNumber2,run);
				delete dataLoopGas2;
				} else{ gasCount = 0;}
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
			} // Loop on Bin
			genLineShape1->Reset("ICESM");
			genLineShape2->Reset("ICESM");
		} // Loop on Repetition
	
		int Tot1 = std::accumulate(v1Tot.begin(), v1Tot.end(), 0);
		int Tot2 = std::accumulate(v2Tot.begin(), v2Tot.end(), 0);
		// SAVE SIMULATED EVENTS 
		ROOT::RDataFrame d1(Tot1 -1);	// PDF1
		ROOT::RDataFrame d2(Tot2 -1); 	// PDF2
		
		TString nameFile1 = TString::Format("LoopDataPdf1_%d.root", l);
		TString nameFile2 = TString::Format("LoopDataPdf2_%d.root", l);
		
		vector<Double_t> rn; vector<Double_t> rn2;
		for(int i = 0; i < v1Type.size(); i++){ rn.push_back(r->Uniform(1));}
		for(int i = 0; i < v2Type.size(); i++){ rn2.push_back(r->Uniform(1));}
		
		int j(0);
		auto FilledFrame1 = FillDataFrame(d1,
						dataPdf1,
						f1,
						v1Type,
						v1Tot,
						j,
						rn,
						RunNumber1,
						cb_shift); // Fill RDataFrame
		std::cout << "Save " <<  folder + nameFile1 << std::endl;
		FilledFrame1.Snapshot("myTree", folder + nameFile1);

		j = 0; // FROM HERE FILL DA DATA
		auto FilledFrame2 = FillDataFrame(d2,
						dataPdf2,
						f2,
						v2Type,
						v2Tot,
						j,
						rn2,
						RunNumber2,
						da_shift); // Fill RDataFrame
		std::cout << "Save " << folder + nameFile2 << std::endl;
		FilledFrame2.Snapshot("myTree", folder + nameFile2);
		//std::cout << "Efficiency: " << Params.Efficiency << " Background: " << Params.CosmicRate << std::endl;
		//std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
		//std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl;
		dataPdf1.Delete(); dataPdf2.Delete();
	} // Loop Trials
	
	// Show the lineshape
	/*
	SetContent(genLineShape1,Nbin1,Cruijff, x_cb_start, x_cb_peak, sigma0_cb, sigma1_cb, k0_cb, k1_cb, Norm_cb);
	SetContent(genLineShape2,Nbin2,Cruijff, x_da_start, x_da_peak, sigma0_da, sigma1_da, k0_da, k1_da, Norm_da);
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
	delete genLineShape1;
	delete genLineShape2;
}
