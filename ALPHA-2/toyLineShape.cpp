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

void toyLineShape(){

int Nbin = 30; 		// Number of Bins
int Ntot = 1000;	// Number of Total Events
double pMix = 0;	// Weight Mix

double pGas = 1 - pMix; // Weight Gas
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

TRandom3 r;

RooRealVar x("x","r [cm]",0.,4.);
//x.setBins(Nbin);

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);

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
RooRealVar Nuw ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk ("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model", "model", RooArgList{gauss_Mix}, RooArgList{Nmix});
RooAddPdf genGas("model1", "model1", RooArgList{Rayleigh}, RooArgList{Nuw});
RooAddPdf genCosmic("model3", "model3", RooArgList{linearFit}, RooArgList{Nbk});

vector<double> f1; 	// Frequencies pdf1
vector<double> v1Nmix;	// Counts Nmix 
vector<double> v1Ngas;	// Counts Ngas
vector<double> v1Tot;

vector<double> f2;	// Frequencies pdf2
vector<double> v2Nmix;  // Counts Nmix
vector<double> v2Ngas;	// Counts Ngas
vector<double> v2Tot;


RooDataSet data("data", "data", RooArgSet(x));

// LOOP, Assign the Nmix and Ngas per frequence
for(int i = 1; i <= Nbin; i++){
	//Pdf1
	Double_t width = histpdf1->GetBinWidth(i);
	Double_t prob = histpdf1->GetBinContent(i);
	prob = prob*width; 			// Probability of the bin
	Nmix.setVal((pMix*Ntot)*prob); 		// Set Nmix
	Nuw.setVal((pGas*Ntot)/Nbin);		// Set Ngas
	int gasCount = 0; int mixCount = 0;
	if(pMix != 0){ // If Nmix different from 0
		RooDataSet *dataLoopMix = genMix.generate(x,Extended()); // Generate data
		data.append(*dataLoopMix); 				 // Append to big dataset
		mixCount = dataLoopMix->sumEntries();
		v1Nmix.push_back(dataLoopMix->sumEntries()); 	// save counts
	}else{
		v1Nmix.push_back(0);
	}
	
	if( pGas != 0){ // If Ngas different from 0
		RooDataSet *dataLoopGas = genGas.generate(x,Extended());
		data.append(*dataLoopGas);
		v1Ngas.push_back(dataLoopGas->sumEntries());
		gasCount = dataLoopGas->sumEntries();
		}
	else{
		v1Ngas.push_back(0);
	}
	f1.push_back(histpdf1->GetBinCenter(i)); 	// save frequency
	v1Tot.push_back(mixCount + gasCount);
	//Pdf2
	width = histpdf2->GetBinWidth(i);
	prob = histpdf2->GetBinContent(i);
	prob = prob*width;
	f2.push_back(histpdf2->GetBinCenter(i));
	v2Nmix.push_back(r.Poisson((pMix*Ntot)*prob));
	v2Tot.push_back(r.Poisson((pMix*Ntot)*prob) + r.Poisson((pGas*Ntot)/Nbin));
}

int trueTot = 0; 
for(int i = 0; i < v1Nmix.size(); i++){ // Total number of events generated from pdf1
	trueTot += v1Nmix[i];
	trueTot += v1Ngas[i];
	if( i == v1Nmix.size()){
	std::cout << "Events generated : " << trueTot << std::endl;
	}
}


// Save the data in RDataFrame
ROOT::RDataFrame d(trueTot-1);
int j(0); 	// Variable for loop
int k(1); 	// Inner Loop, Events belonging to a single frequence
int bin(1);	// Bin number 
TString datafileName = TString::Format("LineShape/ToyShape1f%d.root", static_cast<int>(pMix*100));
d.Define("id", [&j](){ // Id of the events
		return j; 
		})
	.Define("frequence",	// Frequence of the event
	[&bin, &histpdf1](){
		return histpdf1->GetBinCenter(bin); 
		})
	.Define("Type", //
		[&v1Nmix, &v1Ngas, &k, &bin, &Nbin](){
		
		while(v1Nmix[bin-1] + v1Ngas[bin-1] <= 0 && (bin-1) < v1Ngas.size()){ // Check Bin empty
			++bin;
		}
			if(k <= v1Nmix[bin-1]){
				return 0;
			}
			else if( k - v1Nmix[bin-1] > 0){
			std::cout << k <<   " v1Nmix " << v1Nmix[bin-1] << std::endl;
				return 1;
			}
		})
	.Define("radius",	// Generated radius
		[&j, &v1Nmix, &v1Ngas, &data, &bin, &k, &Nbin](){
		if(k >=  v1Nmix[bin-1] + v1Ngas[bin-1]){
			++bin; 	// All counts per freq. are saved, update the bin
			k = 1;	// Set k to 0 for the next frequence inner loop
		}else{
			++k;	// Update inner loop
		}
		std::cout << "Event id: " << j << std::endl; 
		std::cout << "Bin: " << bin <<std::endl;
		std::cout << "Counts per frequence " << v1Nmix[bin-1] + v1Ngas[bin-1] << " and k: " << k << std::endl; 
		
		//data.get(j)->Print("V");
		const RooArgSet &argSet = *(data.get(j));
		++j;		// Update event id
		return static_cast<RooAbsReal&>(argSet["x"]).getVal();
		})
	.Snapshot("pdf1", datafileName);


auto a3 = new TCanvas("a3", "Mix versus Frequencies Normalized");

for(int i = 0; i< f1.size(); i++){
//v1Nmix[i] = v1Nmix[i]/Ntot;
//v2Nmix[i] = v2Nmix[i]/Ntot;
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
