#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "../Headers/toyLineShape.h"

using namespace RooFit;

void toyLineShape(double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5, int NBin = 30, int NTOT = 10000){
	/* Parameters of the Simulation */
int Nbin = NBin;		// Number of Bins
int Ntot = NTOT;		// Number of Total Events
double Ncosmic = (0.492 * Nbin);// Number of Cosmic Events
double pWall_c = Mix_c;		// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Mix_d;		// Weight annihilation on walls for pdf2 (transition d -> a)
double c = C;			// Percentage of division two datasets
	/*				*/

double d = 1 - c; double Nc = Ntot*c; double Nd = Ntot*d;
double pGas_d = 1 - pWall_d; double pGas_c = 1 - pWall_c;

gInterpreter->GenerateDictionary("ToyLine", "../Headers/toyLineShape.h");
TNtuple file_pdf1("pdf1", "pdf1","x:y");
TNtuple file_pdf2("pdf2", "pdf2","x:y");
file_pdf1.ReadFile("lineShape1.csv");
file_pdf2.ReadFile("lineShape2.csv");

vector<double> v1; // Frequency pd1
vector<double> v2; // Counts per frequence pdf1
vector<double> t1; // Frequence pdf2
vector<double> t2; // Counts per frequence pdf2

for (int row = 0; row < file_pdf1.GetEntries(); ++row) {
	file_pdf1.GetEntry(row);
	v1.push_back(file_pdf1.GetArgs()[0]); 	// Extract frequence
	v2.push_back(file_pdf1.GetArgs()[1]);	// Extract Nls
	}

for (int row = 0; row < file_pdf2.GetEntries(); ++row) {
	file_pdf2.GetEntry(row);
	t1.push_back(file_pdf2.GetArgs()[0]); 	// Extract frequence
	t2.push_back(file_pdf2.GetArgs()[1]);	// Extract Nls
	}

Double_t frequence[v1.size()]; Double_t pdf1[v2.size()];
Double_t frequence2[t1.size()]; Double_t pdf2[t2.size()];

std::copy(v1.begin(),v1.end(),frequence); std::copy(v2.begin(),v2.end(),pdf1);
std::copy(t1.begin(),t1.end(),frequence2); std::copy(t2.begin(),t2.end(),pdf2);

// Interpolate the data with spline
TSpline3 *spline1 = new TSpline3("LineShape1", frequence, pdf1, v1.size());
TSpline3 *spline2 = new TSpline3("LineShape2", frequence2, pdf2, t1.size());
TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin,frequence[0], frequence[v1.size() -1]);
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin,frequence2[0], frequence2[t1.size() -1]);

// Set Content Histograms
SetContent(histpdf1,Nbin,spline1);
SetContent(histpdf2,Nbin,spline2);
// Normalize histograms
Double_t factor = 1.;
SetNormalization(histpdf1);
SetNormalization(histpdf2);

RooRealVar x("x","r [cm]",0.,4.);
//x.setBins(Nbin);

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

//MIXING analytic model
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

//PDF UWlosses
//FIXING sigRay to the value from the fit to the model (see analysisMLE.cpp)
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

//PDF Cosmic
//RooRealVar q("q", "q",0.25, 0.25 - 200);
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));

RooRealVar Nmix("Nmix","Nmix",100, -3000, +1000000);
RooRealVar Ngas ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk ("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model", "model", RooArgList{gauss_Mix}, RooArgList{Nmix});
RooAddPdf genGas("model1", "model1", RooArgList{Rayleigh}, RooArgList{Ngas});
RooAddPdf genCosmic("model3", "model3", RooArgList{linearFit}, RooArgList{Nbk});

vector<double> f1; 	// Frequencies pdf1
vector<double> v1Nmix;	// Counts Nmix 
vector<double> v1Ngas;	// Counts Ngas
vector<double> v1Nbk;	// Cosmic Events
vector<double> v1Tot;	// Total Counts
vector<int>    v1Type;	// Type of Event

vector<double> f2;	// Frequencies pdf2
vector<double> v2Nmix;  // Counts Nmix
vector<double> v2Ngas;	// Counts Ngas
vector<double> v2Nbk;	// Cosmic events
vector<double> v2Tot;	// Total Counts
vector<int>    v2Type;	// Type of Event

RooDataSet data("data", "data", RooArgSet(x)); // Dataset to store the events
RooDataSet dati("dati", "dati", RooArgSet(x));

// LOOP, Assign the Nmix and Ngas per frequence
int trueTot1 = 0; int trueTot2 = 0;
for(int i = 1; i <= Nbin; i++){
	//Pdf1
	double prob = ComputeProb(histpdf1,i);	// Probability of the bin
	SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
	int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
	// MIX
	RooDataSet *dataLoopWall = genMix.generate(x,Extended());
	if(dataLoopWall){			// If Nmix different from 0
		data.append(*dataLoopWall);	// Append to big dataset
		SetVectors(dataLoopWall, v1Nmix, v1Type, mixCount, 0);
	}else {v1Nmix.push_back(0);}
	//GAS 
	RooDataSet *dataLoopGas = genGas.generate(x,Extended());
	if(dataLoopGas){			// If Ngas different from 0
		data.append(*dataLoopGas);
		SetVectors(dataLoopGas,v1Ngas,v1Type,gasCount,1);
	}else{ v1Ngas.push_back(0);}
	//COSMIC
	RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
	if(dataCosmic){
		data.append(*dataCosmic);
		SetVectors(dataCosmic,v1Nbk,v1Type,CosmicCount,2);
	}else{v1Nbk.push_back(0);}
	
	f1.push_back(histpdf1->GetBinCenter(i)); 	// save frequency
	v1Tot.push_back(mixCount + gasCount + CosmicCount);
	trueTot1 += mixCount + gasCount + CosmicCount;
	
	//Pdf2
	prob = ComputeProb(histpdf2,i);
	SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
	// MIX
	RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());// Generate data
	if(dataLoopWall){	// Append to big dataset
		dati.append(*dataLoopWall2);		
		SetVectors(dataLoopWall2, v2Nmix, v2Type, mixCount, 0);
	}else{ v2Nmix.push_back(0);}
	//GAS
	RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());
	if(dataLoopGas){	// If Ngas different from 0
		dati.append(*dataLoopGas2);
		SetVectors(dataLoopGas2,v2Ngas,v2Type,gasCount,1);
	}else{ v2Ngas.push_back(0);}
	//COSMIC
	if(dataCosmic){
		dati.append(*dataCosmic);
		SetVectors(dataCosmic,v2Nbk,v2Type,CosmicCount,2);
	}else{v2Nbk.push_back(0);}
	
	f2.push_back(histpdf2->GetBinCenter(i));
	v2Tot.push_back(mixCount + gasCount + CosmicCount);
	trueTot2 += mixCount + gasCount + CosmicCount;
}

//Create the dataframe
ROOT::RDataFrame d1(trueTot1-1); // PDF1
ROOT::RDataFrame d2(trueTot2-1); // PDF2

TString datafileName = TString::Format("ToyShape1_%d_%d_c%d.root", static_cast<int>(pWall_c*100),
static_cast<int>(pWall_d*100),
static_cast<int>(c*100));
TString datafile = TString::Format("ToyShape2_%d_%d_c%d.root", static_cast<int>(pWall_c*100),
static_cast<int>(pWall_d*100),
static_cast<int>(c*100));

FillDataFrame(d1,datafileName, histpdf1, data, v1Type, v1Tot);
FillDataFrame(d2,datafile, histpdf2, dati, v2Type, v2Tot);

auto a3 = new TCanvas("a3", "Mix versus Frequencies Normalized");

for(int i = 0; i< f1.size(); i++){
//v1Nmix[i] = v1Nmix[i]/Nc;
//v2Nmix[i] = v2Nmix[i]/Nc;
}

auto g1 = new TGraph(f1.size(),f1.data(),v1Nmix.data());
auto g2 = new TGraph(f2.size(),f2.data(),v2Nmix.data());
auto pad2 = new TPad("pad2", "pad",0,0,1,1);
pad2->Divide(2,1,0.001,0.001); pad2->Draw();
pad2->cd(1);
g1->SetMarkerStyle(21);
g1->SetTitle("Pdf 1, mix lineshape");
g1->GetYaxis()->SetTitle("Counts");
g1->GetXaxis()->SetTitle("frequencies [MHz]");
g1->Draw();
pad2->cd(2);
g2->SetMarkerStyle(21);
g2->SetTitle("Pdf 2, mix lineshape");
g2->GetYaxis()->SetTitle("Counts");
g2->GetXaxis()->SetTitle("frequencies [MHz]");
g2->Draw();

auto a4 = new TCanvas("a4", "Counts versus Frequencies");
auto g3 = new TGraph(f1.size(),f1.data(),v1Tot.data());
auto g4 = new TGraph(f2.size(),f2.data(),v2Tot.data());
auto pad3 = new TPad("pad2", "pad",0,0,1,1);
pad3->Divide(2,1,0.001,0.001); pad3->Draw();
pad3->cd(1);
g3->SetMarkerStyle(21);
g3->SetTitle("Pdf 1, total count");
g3->GetYaxis()->SetTitle("Counts");
g3->GetXaxis()->SetTitle("frequencies [MHz]");
g3->Draw();
pad3->cd(2);
g4->SetMarkerStyle(21);
g4->SetTitle("Pdf 2, total count");
g4->GetYaxis()->SetTitle("Counts");
g4->GetXaxis()->SetTitle("frequencies [MHz]");
g4->Draw();
}
