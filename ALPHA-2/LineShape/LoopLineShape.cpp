#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "../Headers/toyLineShape.h"

using namespace RooFit;

double model1(double x){
	if(x <= 0){ return 0;}
	else if( x > 0 && x <= 1){ return 2*x;}
	else {return 0;}
}

double model2(double x){
	if(x <= 1420){ return 0;}
	else if( x > 1420 && x <= 1422){ return 2*x;}
	else {return 0;}
}

void LoopLineShape(double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5, int NBin = 30, int NTOT = 5000, int Nloop = 5, bool Save = false, bool MethodSpline = false){
	/* Parameters of the Simulation */
int Nbin = NBin;		// Number of Bins
int Ntot = NTOT;		// Number of Total Events
double Ncosmic = (0.492 * Nbin);// Number of Cosmic Events
double pWall_c = Mix_c;	// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Mix_d;	// Weight annihilation on walls for pdf2 (transition d -> a)
double c = C;			// Percentage of division two datasets
	/*				*/
int Ntrial = Nloop;
double d = 1 - c; double Nc = Ntot*c; double Nd = Ntot*d;
double pGas_d = 1 - pWall_d; double pGas_c = 1 - pWall_c;

gInterpreter->GenerateDictionary("ToyLine","../Headers/toyLineShape.h");

//Histogram to store the distributions
TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin, -1.2096, 1.80480 );
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin, 1419.20962, 1423.86533);

if(MethodSpline){
	SplineMethod(histpdf1,histpdf2, Nbin);
}

if(!MethodSpline){
	SetContent(histpdf1,Nbin,model1);
	SetContent(histpdf2,Nbin,model2);
}

SetNormalization(histpdf1);
SetNormalization(histpdf2);

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

RooRealVar x("x","r [cm]",0.,4.);
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix); //Pdf Annihilation on walls
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay)); //PDF Residual Gas
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x)); //PDF Cosmic

RooRealVar Nmix("Nmix","Nmix",100, -3000, +1000000);
RooRealVar Ngas ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk ("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model", "model", RooArgList{gauss_Mix}, RooArgList{Nmix});
RooAddPdf genGas("model1", "model1", RooArgList{Rayleigh}, RooArgList{Ngas});
RooAddPdf genCosmic("model3", "model3", RooArgList{linearFit}, RooArgList{Nbk});

vector<double> f1, f2;			// Frequencies pdf1 , Frequencies pdf2
vector<double> v1Nmix, v2Nmix;	// Counts Nmix 
vector<double> v1Ngas, v2Ngas;	// Counts Ngas
vector<double> v1Nbk, v2Nbk;	// Cosmic Events
vector<double> v1Tot, v2Tot;	// Total Counts
vector<int>    v1Type,v2Type;	// Type of Event

RooDataSet dataPdf1("data", "data", RooArgSet(x)); // Dataset to store the events
RooDataSet dataPdf2("dati", "dati", RooArgSet(x));

//External Toy Loop
	for(int l = 0; l < Ntrial; l++){
		for(int i = 1; i <= Nbin; i++){ // Inner loop, assign counts to each frequency
			//Pdf1
			double prob = ComputeProb(histpdf1,i);	// Probability of the bin
			SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
			
			if(prob > 0){
			RooDataSet *dataLoopWall = genMix.generate(x,Extended()); // GENERATE THE DATA
			SetVectors(dataPdf1,dataLoopWall,v1Nmix, v1Type, mixCount,   f1,histpdf1->GetBinCenter(i),0); // FILL DATASET
			}else { mixCount = 0; v1Nmix.push_back(0);}
			
			RooDataSet *dataLoopGas = genGas.generate(x, Extended());
			RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
			SetVectors(dataPdf1,dataLoopGas, v1Ngas, v1Type, gasCount,   f1,histpdf1->GetBinCenter(i),1);
			SetVectors(dataPdf1,dataCosmic,  v1Nbk,  v1Type, CosmicCount,f1,histpdf1->GetBinCenter(i),2);
			
			// STORE SOME USEFUL QUANTITIES
			v1Tot.push_back(mixCount + gasCount + CosmicCount);
			
			mixCount = 0; gasCount = 0; CosmicCount = 0;
			// PDF 2
			prob = ComputeProb(histpdf2,i);	// Probability of the bin
			SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			
			if(prob > 0){
			RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
			SetVectors(dataPdf2,dataLoopWall2,v2Nmix, v2Type, mixCount,   f2,histpdf2->GetBinCenter(i),0);
			} else{  mixCount = 0 ; v2Nmix.push_back(0);}
			RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
			// FILL THE DATASET
			
			SetVectors(dataPdf2,dataLoopGas2, v2Ngas, v2Type, gasCount,   f2,histpdf2->GetBinCenter(i),1);
			SetVectors(dataPdf2,dataCosmic,   v2Nbk,  v2Type, CosmicCount,f2,histpdf2->GetBinCenter(i),2);
			// STORE USEFUL QUANTITIES
			v2Tot.push_back(mixCount + gasCount + CosmicCount);
		}
	
	int Tot1 = std::accumulate(v1Tot.begin(), v1Tot.end(), 0);
	int Tot2 = std::accumulate(v2Tot.begin(), v2Tot.end(), 0);
	std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
	std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl;
	
	// DATA ANALYSIS ON SIMULATED EVENTS
	//ROOT::RDataFrame d1(Tot1 -1); // PDF1
	ROOT::RDataFrame d1(Tot1 -1);
	ROOT::RDataFrame d2(Tot2 -1); // PDF2
	int j(0);
	auto FilledFrame1 = FillDataFrame(d1, dataPdf1, f1, v1Type, v1Tot,j); // Fill RDataFrame
	FilledFrame1.Snapshot("myTree", "prova1.root");
	j = 0;
	auto FilledFrame2 = FillDataFrame(d2, dataPdf2, f2, v2Type, v2Tot,j); // Fill RDataFrame
	FilledFrame2.Snapshot("myTree", "prova2.root");
	} // Loop Trials
	
} // End program
