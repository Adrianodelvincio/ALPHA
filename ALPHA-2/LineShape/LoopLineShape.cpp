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

using namespace RooFit;

double linear(double x, double start = 0, double peak = 0,double end = 0){
	//Condizioni raccordo
	double a = (start - peak)/(end - peak);
	double b = (peak - start)/(1 - peak/end);
	if(x <= start){ return 0;}
	else if(x > start && x <= peak){ return (x - start);}
	else if(x > peak && x <= end) { return a*x + b;}
	else { return 0;}
}

double parabola(double x, double start = 0, double peak = 0, double end = 0){
	if(x <= start){ return 0;}
	else if( x > start && x <= peak){ return pow(x - start,2);}
	else if( x > peak && x<= end) {
	double m = -pow(peak - start,2)/(end - peak);
	return m*(x - peak) + pow(peak - start,2);}
	else{return 0;}
}

double Cruijff(double x, double start = 0, double peak){
	double arg = 0;
	

}

void LoopLineShape(	int Nloop = 1,
			TString folder = "linear/",
			TString ConfFile = "ToyConfiguration.txt",					
			bool MethodSpline = false){

gInterpreter->GenerateDictionary("ToyLine","../Headers/toyLineShape.h");
gInterpreter->GenerateDictionary("ToyParser","../Headers/ConfigurationParser.h");
ROOT::EnableImplicitMT(20);
// ReadConfigurationFiles;
ReadConfFile Params(ConfFile);
Params.Print();

//	 Parameters of the Simulation 
int Ntot = Params.Nstack * Params.NHbar * Params.Efficiency;		// Number of Total Events
double FrequencyStep = Params.FrequencyStep;				// Kilo hertz
int SweepStep = Params.SweepStep;					// Number of FrequencyStepS
int Repetition = Params.Repetition;					// Repetition of the single run
double pWall_c = Params.pwall_cb;					// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Params.pwall_ad;					// Weight annihilation on walls for pdf2 (transition d -> a)
double Ncosmic = Params.TimeStep * Params.CosmicRate;			// Number of Cosmic Events
double x_cb_start = Params.x_cb_start;
double x_cb_end = Params.x_cb_end;
double x_cb_peak = Params.x_cb_peak;
double x_da_start = Params.x_da_start;
double x_da_end = Params.x_da_end;
double x_da_peak = Params.x_da_peak;
double RangeDelay = Params.delay;							// Set range delay

	
double startPdf1 = Params.x_cb_start - (FrequencyStep)*5.5;	// Start of frequency sweep c-b
double startPdf2 = Params.x_da_start - (FrequencyStep)*5.5;	// Start of frequency sweep d-a
int Ntrial = Nloop;										// Ntrial
double Nc = Ntot*Params.C; 									// Expected event for lineshape1
double Nd = Ntot*(1 - Params.C); 								// Expected event for lineshape2
double pGas_d = 1 - pWall_d; 									// Percentage of annihilation on residual gas
double pGas_c = 1 - pWall_c; 									// Percentage of annihilation on residual gas

int Nbin1 = SweepStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP
int Nbin2 = SweepStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP

//if(MethodSpline){ SplineMethod(histpdf1,histpdf2, Nbin);}


RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

RooRealVar x("x","r [cm]",0.,4.);
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
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

TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin1, startPdf1, startPdf1 + Nbin1*(FrequencyStep));
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin2, startPdf2, startPdf2 + Nbin2*(FrequencyStep));

//External Toy Loop
TRandom3 *r = new TRandom3();

	for(int l = 0; l < Ntrial; l++){				//LOOP ON TRIALS
		RooDataSet dataPdf1("data", "data", RooArgSet(x));	// Dataset to store the events
		RooDataSet dataPdf2("dati", "dati", RooArgSet(x));	// Dataset to store the events
		vector<double>	f1, f2;					// Frequencies pdf1 , Frequencies pdf2
		vector<double>	v1Tot, v2Tot;				// Total Counts
		vector<int>	v1Type,v2Type;				// Type of Event
		vector<int>	RunNumber1, RunNumber2;			// Run Number (from 0 to Repetition)
		vector<double>  freqDelay;
		for(int run = 0; run < Repetition; run++){		// LOOP ON REPETITION
			//DEFINITION OF LINESHAPE
			// SET DELAY FOR THE ONSET
			double delay = r->Uniform(-RangeDelay,RangeDelay);
			//double delay = RangeDelay;
			freqDelay.push_back(delay);
			double start1 = x_cb_start + delay;
			double start2 = x_da_start + delay;
			// SET CONTENT OF THE HISTOGRAMS AND NORMALIZE
			SetContent(histpdf1,Nbin1,parabola, start1,x_cb_peak ,x_cb_end);
			SetContent(histpdf2,Nbin2,parabola, start2,x_da_peak ,x_da_end);
			SetNormalization(histpdf1);
			SetNormalization(histpdf2);
			
			for(int bin = 1; bin <= SweepStep; bin++){ 	//LOOP ON BINS
				//PDF1
				double prob = ComputeProb(histpdf1,bin);// Probability of the bin
				SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin1,Ncosmic, &Nwall,&Ngas,&Nbk);
				int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
				double frequence = histpdf1->GetBinCenter(bin);
				
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
				
				// PDF 2
				mixCount = 0; gasCount = 0; CosmicCount = 0;
				prob = ComputeProb(histpdf2,bin);	// Probability of the bin
				frequence = histpdf2->GetBinCenter(bin);
				SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin2,Ncosmic, &Nwall,&Ngas,&Nbk);
				
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
			} // Loop on Bin
			if(l != Ntrial-1 && run != Repetition -1){
			histpdf1->Reset("ICESM");
			histpdf2->Reset("ICESM");}
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
						freqDelay); // Fill RDataFrame
		std::cout << "Save " <<  folder + nameFile1 << std::endl;
		FilledFrame1.Snapshot("myTree", folder + nameFile1);

		j = 0; // FROM HERE APPLY THE ALGORITHM TO PDF2
		auto FilledFrame2 = FillDataFrame(d2,
						dataPdf2,
						f2,
						v2Type,
						v2Tot,
						j,
						rn2,
						RunNumber2,
						freqDelay); // Fill RDataFrame
		std::cout << "Save " << folder + nameFile2 << std::endl;
		FilledFrame2.Snapshot("myTree", folder + nameFile2);		
		std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
		std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl;
		dataPdf1.Delete(); dataPdf2.Delete();
	} // Loop Trials
	
	auto b = new TCanvas("b1", "Spectral lines");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	gStyle->SetOptStat(0);
	histpdf1->GetXaxis()->SetTitle("frequency [kHz]");
	histpdf1->SetMarkerStyle(21);
	histpdf1->SetMarkerColor(2);
	histpdf1->SetLineColor(4);
	histpdf1->Draw();
	pad->cd(2);
	histpdf2->GetXaxis()->SetTitle("frequency [kHz]");
	histpdf2->SetMarkerStyle(21);
	histpdf2->SetMarkerColor(2);
	histpdf2->SetLineColor(4);
	histpdf2->Draw();
} // End program

