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

double linear(double x, double start = 0, double peak = 0,double end = 0){
	//Condizioni raccordo
	double a = (start - peak)/(end - peak);
	double b = (peak - start)/(1 - peak/end);
	if(x <= start){ return 0;}
	else if(x > start && x <= peak){ return (x - start);}
	else if(x > peak && x <= end) { return a*x + b;}
	else { return 0;}
}

double parabola1(double x, double xmin = 0, double xmax = 0){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return 3*(x*x);}
	else {return 0;}
}

void LoopLineShape(TString folder = "Linear/",
				int Nloop = 1,
				double mix_cb = 1,
				double mix_ad = 1,
				double C = 0.5,
				int Nstack = 20,
				int NHbar = 14,
				int Repetition = 5,
				double FrequencyStep = 5,
				int timeStep = 8,
				int SweepStep = 24,
				double x_cb_start = 0,
				double x_cb_peak = 20,
				double x_cb_end = 40,
				double x_da_start = 1420000,
				double x_da_peak = 1420020,	
				double x_da_end = 1420040,
				double CosmicRate = 1786./35000.,
				double Efficiency = 1,
				bool MethodSpline = false){

	/* Parameters of the Simulation */
int Ntot = Nstack * NHbar * Efficiency;	// Number of Total Events
double pWall_c = mix_cb;				// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = mix_ad;				// Weight annihilation on walls for pdf2 (transition d -> a)
double c = C;							// Percentage of division two datasets
FrequencyStep = FrequencyStep;			// Kilo hertz
	/*				*/

int Ntrial = Nloop;									// Ntrial
double d = 1 - c; 									// Percentage of event for lineshape1
double Nc = Ntot*c; 								// Expected event for lineshape1
double Nd = Ntot*d; 								// Expected event for lineshape2
double pGas_d = 1 - pWall_d; 						// Percentage of annihilation on residual gas
double pGas_c = 1 - pWall_c; 						// Percentage of annihilation on residual gas
double Ncosmic = SweepStep * timeStep * CosmicRate;	// Number of Cosmic Events
double startPdf1 = x_cb_start -(FrequencyStep)*5;	// Start of frequency sweep c-b
double startPdf2 = x_da_start -(FrequencyStep)*5;	// Start of frequency sweep d-a


gInterpreter->GenerateDictionary("ToyLine","../Headers/toyLineShape.h");
ROOT::EnableImplicitMT();

int Nbin1 = SweepStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP
int Nbin2 = SweepStep;	// FIX NUMBER BIN EQUAL TO SWEEPSTEP
//Histogram to store the distributions
//TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin1, startPdf1, startPdf1 + Nbin1*(FrequencyStep));
//TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin2, startPdf2, startPdf2 + Nbin2*(FrequencyStep));

/*if(MethodSpline){ // Set the lineshape with the Spline method
	SplineMethod(histpdf1,histpdf2, Nbin);
}

if(!MethodSpline){ // DEFINE CUSTOM LINESHAPES
	SetContent(histpdf1,Nbin1,linear, x_cb_start,x_cb_peak ,x_cb_end);
	SetContent(histpdf2,Nbin2,linear, x_da_start,x_da_peak ,x_da_end);
}

SetNormalization(histpdf1);
SetNormalization(histpdf2);
*/
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

	for(int l = 0; l < Ntrial; l++){					//LOOP ON TRIALS
		RooDataSet dataPdf1("data", "data", RooArgSet(x));	// Dataset to store the events
		RooDataSet dataPdf2("dati", "dati", RooArgSet(x));	// Dataset to store the events
		vector<double>	f1, f2;								// Frequencies pdf1 , Frequencies pdf2
		vector<double>	v1Tot, v2Tot;						// Total Counts
		vector<int>		v1Type,v2Type;						// Type of Event
		vector<int>		RunNumber1, RunNumber2;				// Run Number (from 0 to Repetition)
		vector<double>  freqDelay;
		for(int run = 0; run < Repetition; run++){		// LOOP ON REPETITION
			//DEFINITION OF LINESHAPE
			// SET DELAY FOR THE ONSET
			double delay = r->Uniform(-FrequencyStep/2, +FrequencyStep/2);
			freqDelay.push_back(delay);
			double start1 = x_cb_start + delay;
			double start2 = x_da_start + delay;
			// SET CONTENT OF THE HISTOGRAMS AND NORMALIZE
			SetContent(histpdf1,Nbin1,linear, start1,x_cb_peak ,x_cb_end);
			SetContent(histpdf2,Nbin2,linear, start2,x_da_peak ,x_da_end);
			SetNormalization(histpdf1);
			SetNormalization(histpdf2);
			
			for(int bin = 1; bin <= SweepStep; bin++){ 	//LOOP ON BINS
				//PDF1
				double prob = ComputeProb(histpdf1,bin);	// Probability of the bin
				SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin1,Ncosmic/Nbin1, &Nwall,&Ngas,&Nbk);
				int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
				int frequence = histpdf1->GetBinCenter(bin);
				if(prob > 0){
				RooDataSet *dataLoopWall = genMix.generate(x,Extended()); // GENERATE THE DATA
				SetVectors(dataPdf1,dataLoopWall,
							v1Type, 0,
							mixCount,
							f1,frequence,
							RunNumber1,run); // FILL DATASET
				}else { mixCount = 0;}
				
				if(Ngas.getVal() != 0){
				RooDataSet *dataLoopGas = genGas.generate(x, Extended());
				SetVectors(dataPdf1,dataLoopGas,
							v1Type, 1,
							gasCount,
							f1,frequence,
							RunNumber1,run);
				}else {gasCount = 0;}
				
				RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
				SetVectors(dataPdf1,dataCosmic,
							v1Type, 2,
							CosmicCount,
							f1,frequence,
							RunNumber1,run);
				
				// STORE SOME USEFUL QUANTITIES
				v1Tot.push_back(mixCount + gasCount + CosmicCount);
				
				// PDF 2
				mixCount = 0; gasCount = 0; CosmicCount = 0;
				prob = ComputeProb(histpdf2,bin);	// Probability of the bin
				frequence = histpdf2->GetBinCenter(bin);
				SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin2,Ncosmic/Nbin2, &Nwall,&Ngas,&Nbk);
				
				if(prob > 0){
				RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
				SetVectors(dataPdf2,dataLoopWall2,
							v2Type,0,
							mixCount,
							f2,frequence,
							RunNumber2,run);
				} else{  mixCount = 0;}
				if(Ngas.getVal() != 0){
				RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
				SetVectors(dataPdf2,dataLoopGas2,
							v2Type,1,
							gasCount,
							f2,frequence,
							RunNumber2,run);
				} else{ gasCount = 0;}
				RooDataSet *dataCosmic2 = genCosmic.generate(x, Extended());			
				SetVectors(dataPdf2,dataCosmic2,
							v2Type,2,
							CosmicCount,
							f2,frequence,
							RunNumber2, run);
				// STORE USEFUL QUANTITIES
				v2Tot.push_back(mixCount + gasCount + CosmicCount);
			} // Loop on Bin
			histpdf1->Reset("ICESM");
			histpdf2->Reset("ICESM");
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

	} // Loop Trials
} // End program

