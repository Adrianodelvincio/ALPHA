#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "../Headers/toyLineShape.h"

using namespace RooFit;

void LoadLineShapeData(TNtuple*,TNtuple*, TString file1 = "lineShape1.csv" , TString file2 = "lineShape2.csv"); 

void LoopLineShape(double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5, int NBin = 30, int NTOT = 5000, int Nloop = 5, bool Save = false){
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

TNtuple *file1;
TNtuple *file2;

LoadLineShapeData(file1,file2);

TNtuple file_pdf1("pdf1", "pdf1","x:y");
TNtuple file_pdf2("pdf2", "pdf2","x:y");
file_pdf1.ReadFile("lineShape1.csv");
file_pdf2.ReadFile("lineShape2.csv");

vector<double> v1,v2,t1,t2; // Frequency pd1, Counts per frequence pdf1, Frequence pdf2,Counts per frequence pdf2

ConvertTNtutpla(file_pdf1,v1,v2);
ConvertTNtutpla(file_pdf2,t1,t2);
Double_t frequence[v1.size()]; Double_t pdf1[v2.size()];
Double_t frequence2[t1.size()]; Double_t pdf2[t2.size()];
std::copy(v1.begin(),v1.end(),frequence); std::copy(v2.begin(),v2.end(),pdf1);
std::copy(t1.begin(),t1.end(),frequence2); std::copy(t2.begin(),t2.end(),pdf2);

// Interpolate the data with spline
TSpline3 *spline1 = new TSpline3("LineShape1", frequence, pdf1, v1.size());
TSpline3 *spline2 = new TSpline3("LineShape2", frequence2, pdf2, t1.size());
TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin,frequence[0], frequence[v1.size() -1]);
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin,frequence2[0], frequence2[t1.size() -1]);

// Set Content Histograms, Normalize histograms
SetContent(histpdf1,Nbin,spline1);
SetContent(histpdf2,Nbin,spline2);
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
		
		// LOOP, Assign the Nmix and Ngas per frequence
		int Tot1 = 0; int Tot2 = 0;
		for(int i = 1; i <= Nbin; i++){
			//Pdf1
			double prob = ComputeProb(histpdf1,i);	// Probability of the bin
			SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
			
			// GENERATE THE DATA
			RooDataSet *dataLoopWall = genMix.generate(x,Extended());
			RooDataSet *dataLoopGas = genGas.generate(x, Extended());
			RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
			// FILL DATASET
			SetVectors(dataPdf1,dataLoopWall,v1Nmix, v1Type, mixCount,   f1,histpdf1->GetBinCenter(i),0);
			SetVectors(dataPdf1,dataLoopGas, v1Ngas, v1Type, gasCount,   f1,histpdf1->GetBinCenter(i),1);
			SetVectors(dataPdf1,dataCosmic,  v1Nbk,  v1Type, CosmicCount,f1,histpdf1->GetBinCenter(i),2);
			// STORE SOME USEFUL QUANTITIES
			v1Tot.push_back(mixCount + gasCount + CosmicCount);
			Tot1 += mixCount + gasCount + CosmicCount;
			mixCount = 0; gasCount = 0; CosmicCount = 0;
			
			// PDF 2
			prob = ComputeProb(histpdf2,i);
			SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			
			// GENERATE THE DATA
			RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
			RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
			// FILL THE DATASET
			SetVectors(dataPdf2,dataLoopWall,v2Nmix, v2Type, mixCount,   f2,histpdf2->GetBinCenter(i),0);
			SetVectors(dataPdf2,dataLoopGas, v2Ngas, v2Type, gasCount,   f2,histpdf2->GetBinCenter(i),1);
			SetVectors(dataPdf2,dataCosmic,  v2Nbk,  v2Type, CosmicCount,f2,histpdf2->GetBinCenter(i),2);
			// STORE USEFUL QUANTITIES
			v2Tot.push_back(mixCount + gasCount + CosmicCount);
			Tot2 += mixCount + gasCount + CosmicCount;
		}
	
	std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
	std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl;
	
	// DATA ANALYSIS ON SIMULATED EVENTS
	//ROOT::RDataFrame d1(Tot1 -1); // PDF1
	ROOT::RDataFrame d1(Tot1 -1);
	ROOT::RDataFrame d2(Tot2 -1); // PDF2
	int j(0);
	auto FilledFrame1 = FillDataFrame(d1, dataPdf1 ,f1,v2Type,v2Tot,j); // Fill RDataFrame
	FilledFrame1.Snapshot("myTree", "prova1.root");
	j = 0;
	auto FilledFrame2 = FillDataFrame(d2, dataPdf2,f2 , v2Type, v2Tot,j); // Fill RDataFrame
	FilledFrame2.Snapshot("myTree", "prova2.root");
	} // Loop Trials
	
} // End program


void LoadLineShapeData(TNtuple *file_pdf1, TNtuple *file_pdf2,TString file1 = "lineShape1.csv" , TString file2 = "lineShape2.csv" ){
	file_pdf1 = new TNtuple("pdf1", "pdf1","x:y");
	file_pdf2 = new TNtuple("pdf2", "pdf2","x:y");
	file_pdf1->ReadFile(file1);
	file_pdf2->ReadFile(file2);
}

