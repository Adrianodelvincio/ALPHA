#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "Headers/toyLineShape.h"
#include "TRandom.h"
using namespace RooFit;

void toyLineShape(int Flag = 0, double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5){
/////
int Nbin = 30; 		// Number of Bins
int Ntot = 10000;	// Number of Total Events
int Ncosmic = static_cast<int>(0.492 * Nbin);	// Number of Cosmic Events
//int Ncosmic = 100;	// Number of Cosmic Events
double pMix_c = 0.5;	// Weight MIx pdf1
double pMix_d = 0.5;	// Weight Mix pdf2
double c = 0.5;		// Percentage of division two datasets
/////

if(Flag){ // Pass from command line
pMix_c = Mix_c;
pMix_d = Mix_d;
c = C;
}

double d = 1 - c;
Ntot = Ntot - Ncosmic;
double Nc = Ntot*c;
double Nd = Ntot*d;
double pGas_d = 1 - pMix_d; // Weight Gas
double pGas_c = 1 - pMix_c;
gInterpreter->GenerateDictionary("ToyLine", "Headers/toyLineShape.h");
TNtuple file_pdf1("pdf1", "pdf1","x:y");
TNtuple file_pdf2("pdf2", "pdf2","x:y");
file_pdf1.ReadFile("LineShape/lineShape1.csv");
file_pdf2.ReadFile("LineShape/lineShape2.csv");

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
histpdf1->Scale(factor/histpdf1->GetEntries());
histpdf2->Scale(factor/histpdf2->GetEntries());

//Visualize spline
TCanvas *a1 = new TCanvas("a1","Spiline interpolation");
auto pad = new TPad("pad", "pad",0,0,1,1);
pad->Divide(2,1,0.,0.); pad->Draw();
pad->cd(1);
file_pdf1.SetMarkerStyle(21);
file_pdf1.Draw("y:x");
spline1->SetLineColor(kRed);
spline1->SetLineWidth(3);
spline1->Draw("same");
pad->cd(2);
file_pdf2.SetMarkerStyle(21);
file_pdf2.Draw("y:x");
spline2->SetLineColor(kRed);
spline2->SetLineWidth(3);
spline2->Draw("same");


TCanvas *a2 = new TCanvas("a2", "Extracted histogram");
auto pad1 = new TPad("pad1", "pad",0,0,1,1);
pad1->Divide(2,1,0.,0.); pad1->Draw();
pad1->cd(1);
histpdf1->Draw();
pad1->cd(2);
histpdf2->Draw();


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

TRandom3 r;
RooDataSet data("data", "data", RooArgSet(x)); // Dataset to store the events
RooDataSet dati("dati", "dati", RooArgSet(x));

// LOOP, Assign the Nmix and Ngas per frequence
for(int i = 1; i <= Nbin; i++){
	//Pdf1
	double prob = ComputeProb(histpdf1,i);	// Probability of the bin
	SetCoefficients((pMix_c*Nc)*prob,(pGas_c*Nc)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
	int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
	// MIX
	RooDataSet *dataLoopMix = genMix.generate(x,Extended());
	if(dataLoopMix){			// If Nmix different from 0
		data.append(*dataLoopMix);	// Append to big dataset
		SetVectors(dataLoopMix, v1Nmix, v1Type, mixCount, 0);
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
	
	//Pdf2
	prob = ComputeProb(histpdf2,i);
	SetCoefficients((pMix_d*Nd)*prob,(pGas_d*Nd)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
	// MIX
	RooDataSet *dataLoopMix2 = genMix.generate(x,Extended());// Generate data
	if(dataLoopMix){	// Append to big dataset
		dati.append(*dataLoopMix2);		
		SetVectors(dataLoopMix2, v2Nmix, v2Type, mixCount, 0);
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
}

int trueTot = 0; 
for(int i = 0; i < v1Nmix.size(); i++){ // Total number of events generated from pdf1
	trueTot += v1Nmix[i];
	trueTot += v1Ngas[i];
	trueTot += v1Nbk[i];
	if( i == v1Nmix.size()){
	std::cout << "Events generated : " << trueTot << std::endl;
	}
}

ROOT::RDataFrame d1(trueTot-1); // PDF1

trueTot = 0;
for(int i = 0; i < v2Nmix.size(); i++){ // Total number of events generated from pdf1
	trueTot += v2Nmix[i];
	trueTot += v2Ngas[i];
	trueTot += v2Nbk[i];
	if( i == v2Nmix.size()){
	std::cout << "Events generated : " << trueTot << std::endl;
	}
}
ROOT::RDataFrame d2(trueTot -1); // PDF2

// 
int j(0); 	// Variable for loop
int k(1); 	// Inner Loop, Events belonging to a single frequence
int bin(1);	// Bin number 
TString datafileName = TString::Format("LineShape/ToyShape1_%d_%d_c%d.root", static_cast<int>(pMix_c*100),
static_cast<int>(pMix_d*100),
static_cast<int>(c*100));

d1.Define("id", [&j](){		// Id of the events
		return j; 
		})
	.Define("frequence",	// Frequence of the event
	[&bin, &histpdf1](){
		return histpdf1->GetBinCenter(bin); 
		})
	.Define("Type", //
		[&v1Type, &j](){
		return v1Type[j];
		})
	.Define("radius",	// Generated radius
		[&j,&v1Tot, &data, &bin, &k](){
		if(k >=  v1Tot[bin-1]){
			++bin; 	// All counts per freq. are saved, update the bin
			k = 1;	// Set k to 0 for the next frequence inner loop
		}else{
			++k;	// Update inner loop
		}
		std::cout << "Event id: " << j << std::endl; 
		std::cout << "Bin: " << bin-1 << " Counts per freq: " << v1Tot[bin -1] << " k event: " << k; 
		data.get(j)->Print("V");
		const RooArgSet &argSet = *(data.get(j));
		++j;		// Update event id
		return static_cast<RooAbsReal&>(argSet["x"]).getVal();
		})
	.Snapshot("myTree", datafileName);

TString datafile = TString::Format("LineShape/ToyShape2_%d_%d_c%d.root", static_cast<int>(pMix_c*100),
static_cast<int>(pMix_d*100),
static_cast<int>(c*100));

j = 0; k = 1; bin = 1;
d2.Define("id", [&j](){		// Id of the events
		return j; 
		})
	.Define("frequence",	// Frequence of the event
	[&bin, &histpdf2](){
		return histpdf2->GetBinCenter(bin); 
		})
	.Define("Type", //
		[&v2Type, &j](){
		return v2Type[j];
		})
	.Define("radius",	// Generated radius
		[&j,&v2Tot, &dati, &bin, &k](){
		if(k >=  v2Tot[bin-1]){
			++bin; 	// All counts per freq. are saved, update the bin
			k = 1;	// Set k to 0 for the next frequence inner loop
		}else{
			++k;	// Update inner loop
		}
		std::cout << "Event id: " << j << std::endl; 
		std::cout << "Bin: " << bin-1 << " Counts per freq: " << v2Tot[bin -1] << " k event: " << k; 
		dati.get(j)->Print("V");
		const RooArgSet &argSet = *(dati.get(j));
		++j;		// Update event id
		return static_cast<RooAbsReal&>(argSet["x"]).getVal();
		})
	.Snapshot("myTree", datafile);


auto a3 = new TCanvas("a3", "Mix versus Frequencies Normalized");

for(int i = 0; i< f1.size(); i++){
//v1Nmix[i] = v1Nmix[i]/Nc;
//v2Nmix[i] = v2Nmix[i]/Nc;
}

auto g1 = new TGraph(f1.size(),f1.data(),v1Nmix.data());
auto g2 = new TGraph(f2.size(),f2.data(),v2Nmix.data());
auto pad2 = new TPad("pad2", "pad",0,0,1,1);
pad2->Divide(2,1,0.,0.); pad2->Draw();
pad2->cd(1);
g1->SetMarkerStyle(21);
g1->SetTitle("MIX");
g1->GetYaxis()->SetTitle("Counts");
g1->GetXaxis()->SetTitle("frequencies [#MHz]");
g1->Draw();
pad2->cd(2);
g2->SetMarkerStyle(21);
g2->SetTitle("MIX");
g2->GetYaxis()->SetTitle("Counts");
g2->GetXaxis()->SetTitle("frequencies [#MHz]");
g2->Draw();

auto a4 = new TCanvas("a4", "Counts versus Frequencies");
auto g3 = new TGraph(f1.size(),f1.data(),v1Tot.data());
auto g4 = new TGraph(f2.size(),f2.data(),v2Tot.data());
auto pad3 = new TPad("pad2", "pad",0,0,1,1);
pad3->Divide(2,1,0.,0.); pad3->Draw();
pad3->cd(1);
g3->SetMarkerStyle(21);
g3->SetTitle("MIX");
g3->GetYaxis()->SetTitle("Counts");
g3->GetXaxis()->SetTitle("frequencies [#MHz]");
g3->Draw();
pad3->cd(2);
g4->SetMarkerStyle(21);
g4->SetTitle("MIX");
g4->GetYaxis()->SetTitle("Counts");
g4->GetXaxis()->SetTitle("frequencies [#MHz]");
g4->Draw();
}
