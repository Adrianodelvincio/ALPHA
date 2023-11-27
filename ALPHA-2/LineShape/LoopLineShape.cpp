#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
//#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include <TRandom3.h>
#include <math.h> 
#include "../Headers/toyLineShape.h"

using namespace RooFit;

double modello1(double x, double xmin = 0, double xmax = 40){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return 2*x;}
	else {return 0;}
}

double modello2(double x, double xmin = 100, double xmax = 140){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return 2*(x - xmin);}
	else {return 0;}
}

double parabola1(double x, double xmin = 0, double xmax = 40){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return 3*(x*x);}
	else {return 0;}
}

double parabola2(double x, double xmin = 100, double xmax = 140){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return (3/2 )*(x - xmin)*(x - xmin);}
	else {return 0;}
}

void LoopLineShape(int Nloop = 400,
				double Mix_c = 1,
				double Mix_d = 1,
				double C = 0.5,
				int Nstack = 8,
				int NHbar = 14,
				int Repetition = 5,
				double FrequencyStep = 5,
				int timeStep = 8,
				int SweepStep = 24,
				double x0 = 0,
				double y0 = 40,
				double x1 = 120,
				double y1 = 160,
				bool MethodSpline = false){

	/* Parameters of the Simulation */
int Ntot = Nstack * NHbar * Repetition;	// Number of Total Events
double pWall_c = Mix_c;				// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Mix_d;				// Weight annihilation on walls for pdf2 (transition d -> a)
double c = C;						// Percentage of division two datasets
FrequencyStep = FrequencyStep;		// in Kilo hertz
	/*				*/

int Ntrial = Nloop;	// Ntrial
double d = 1 - c; 	// Percentage of event for lineshape1
double Nc = Ntot*c; // Expected event for lineshape1
double Nd = Ntot*d; // Expected event for lineshape2
double pGas_d = 1 - pWall_d; // Percentage of annihilation on residual gas
double pGas_c = 1 - pWall_c; // Percentage of annihilation on residual gas
double Ncosmic = SweepStep*Repetition * timeStep * (1786./35000.);	// Number of Cosmic Events
double startPdf1 = -40;	// Range of lineshape1 
double startPdf2 = 100;	// Range of lineshape2


gInterpreter->GenerateDictionary("ToyLine","../Headers/toyLineShape.h");
ROOT::EnableImplicitMT();

int Nbin1 = SweepStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP
int Nbin2 = SweepStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP
//Histogram to store the distributions
TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin1, startPdf1, startPdf1 + Nbin1*(FrequencyStep));
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin2, startPdf2, startPdf2 + Nbin2*(FrequencyStep));

/*if(MethodSpline){ // Set the lineshape with the Spline method
	SplineMethod(histpdf1,histpdf2, Nbin);
}*/

if(!MethodSpline){ // DEFINE CUSTOM LINESHAPES
	SetContent(histpdf1,Nbin1,parabola1, x0, y0);
	SetContent(histpdf2,Nbin2,parabola2, x1, y1);
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

RooRealVar Nwall("Nwall","Nwall",100, -3000, +1000000);
RooRealVar Ngas ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk 	("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model","model", RooArgList{gauss_Mix}, RooArgList{Nwall});
RooAddPdf genGas("model1","model1", RooArgList{Rayleigh}, RooArgList{Ngas});
RooAddPdf genCosmic("model2", "model2", RooArgList{linearFit}, RooArgList{Nbk});

//External Toy Loop
TRandom3 *r = new TRandom3();
	for(int l = 0; l < Ntrial; l++){
		RooDataSet dataPdf1("data", "data", RooArgSet(x));	// Dataset to store the events
		RooDataSet dataPdf2("dati", "dati", RooArgSet(x));	// Dataset to store the events
		vector<double> f1, f2;			// Frequencies pdf1 , Frequencies pdf2
		vector<double> v1Tot, v2Tot;	// Total Counts
		vector<int>    v1Type,v2Type;	// Type of Event
		
		for(int i = 1; i <= SweepStep; i++){ // Inner loop, assign counts to each frequency
			//Pdf1
			double prob = ComputeProb(histpdf1,i);	// Probability of the bin
			SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin1,Ncosmic/Nbin1, &Nwall,&Ngas,&Nbk);
			int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
			
			if(prob > 0){
			RooDataSet *dataLoopWall = genMix.generate(x,Extended()); // GENERATE THE DATA
			SetVectors(dataPdf1,dataLoopWall, v1Type, mixCount,  f1,histpdf1->GetBinCenter(i),0); // FILL DATASET
			}else { mixCount = 0;}
			
			if(Ngas.getVal() != 0){
			RooDataSet *dataLoopGas = genGas.generate(x, Extended());
			SetVectors(dataPdf1,dataLoopGas, v1Type, gasCount,   f1,histpdf1->GetBinCenter(i),1);
			}else {gasCount = 0;}
			
			RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
			SetVectors(dataPdf1,dataCosmic,  v1Type, CosmicCount,f1,histpdf1->GetBinCenter(i),2);
			
			// STORE SOME USEFUL QUANTITIES
			//std::cout << "conteggi residual gas: " << gasCount << std::endl;
			v1Tot.push_back(mixCount + gasCount + CosmicCount);
			
			// PDF 2
			mixCount = 0; gasCount = 0; CosmicCount = 0;
			prob = ComputeProb(histpdf2,i);	// Probability of the bin
			SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin2,Ncosmic/Nbin2, &Nwall,&Ngas,&Nbk);
			
			if(prob > 0){
			RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
			SetVectors(dataPdf2,dataLoopWall2, v2Type, mixCount, f2,histpdf2->GetBinCenter(i),0);
			} else{  mixCount = 0;}
			if(Ngas.getVal() != 0){
			RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
			SetVectors(dataPdf2,dataLoopGas2, v2Type, gasCount,   f2,histpdf2->GetBinCenter(i),1);
			} else{ gasCount = 0;}
			RooDataSet *dataCosmic2 = genCosmic.generate(x, Extended());			
			SetVectors(dataPdf2,dataCosmic2,   v2Type, CosmicCount,f2,histpdf2->GetBinCenter(i),2);
			
			// STORE USEFUL QUANTITIES
			//std::cout << "conteggi residual gas: " << gasCount << std::endl;
			v2Tot.push_back(mixCount + gasCount + CosmicCount);
		}
	
		int Tot1 = std::accumulate(v1Tot.begin(), v1Tot.end(), 0);
		int Tot2 = std::accumulate(v2Tot.begin(), v2Tot.end(), 0);
		// SAVE SIMULATED EVENTS 
		ROOT::RDataFrame d1(Tot1 -1);	// PDF1
		ROOT::RDataFrame d2(Tot2 -1); 	// PDF2
		
		TString nameFile1 = TString::Format("Quadratic/LoopDataPdf1_%d.root", l);
		TString nameFile2 = TString::Format("Quadratic/LoopDataPdf2_%d.root", l);
		
		
		vector<Double_t> rn; vector<Double_t> rn2;
		for(int i = 0; i < v1Type.size(); i++){ rn.push_back(r->Uniform(1));}
		for(int i = 0; i < v2Type.size(); i++){ rn2.push_back(r->Uniform(1));}
		
		int j(0);
		auto FilledFrame1 = FillDataFrame(d1, dataPdf1, f1, v1Type, v1Tot,j,rn); // Fill RDataFrame
		std::cout << "Save " <<  nameFile1 << std::endl;
		FilledFrame1.Snapshot("myTree", nameFile1);

		j = 0; // FROM HERE APPLY THE ALGORITHM TO PDF2
		auto FilledFrame2 = FillDataFrame(d2, dataPdf2, f2, v2Type, v2Tot,j,rn2); // Fill RDataFrame
		std::cout << "Save " << nameFile2 << std::endl;
		FilledFrame2.Snapshot("myTree", nameFile2);

		std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
		std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl; 

	} // Loop Trials
} // End program

