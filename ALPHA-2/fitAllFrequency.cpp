#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
#include "TGraphErrors.h"
#include <string>
#include "Headers/fitAllFrequency.h"
#include <fstream>
using namespace RooFit;

void fitAllFrequency(){
// DATA TO FIT
TString cartella = TString::Format("DataSetROOT/");
TString firstpart = TString::Format("r68481_f");
TString formato = TString::Format(".root");
vector<int> vi;
vector<double> frequency;
vector<int> Counts;
vector<double> resGas;
vector<double> errorx;
vector<double> errGas;
TCanvas *c[8];

//Generate dictionary
gInterpreter->GenerateDictionary("fitAllFrequency", "Headers/fitAllFrequency.h");

//Fit distributions
RooRealVar x("x", "r [cm]", 0, 4);
x.setBins(30);

//RAYLEIGH
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay)))", RooArgSet(x,sigRay));
//COSMIC
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));
//MIXING
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85, 0.85,0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

//COMPLETE MODEL
RooRealVar Nmix("Nfit_{mix}","Nfit_{mix}",0.,200.);
RooRealVar Ngas("Nfit_{gas}","Nfit_{gas}", -100,100);
RooRealVar Ncosmic("Nfit_{cosmic}", "Nfit_{cosmic}",10.2,10.2);
RooAddPdf model_analytic("model", "model", RooArgList{gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix,Ngas,Ncosmic});
RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

int i = 1; // i frequency
while(1){
	
	TString frequence = TString::Format("%d", i); 
	frequency.push_back(i); vi.push_back(i); i += 1;
	std::cout << cartella + firstpart + frequence + formato << std::endl;
	if(gSystem->AccessPathName(cartella + firstpart + frequence + formato)){
		break;
	}
	if(CheckEmpty(cartella + firstpart + frequence + formato)){
		continue;
	}
	//ROOTDATAFRAME MACHINERY
	ROOT::RDataFrame rdf("myTree", cartella + firstpart + frequence + formato);
	auto rdf2 = rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
	auto Hist = rdf2.Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius");
	//Save the entries
	Counts.push_back(Hist->GetEntries());
	resGas.push_back(static_cast<int> (Ngas.getVal()));
	errGas.push_back(static_cast<int> (Ngas.errorVar()->getVal()));
	errorx.push_back(0.);
	RooDataHist hh("dh3", "dh3", x, Import(*Hist));
	RooPlot *analyticframe = x.frame(Title(firstpart +  frequence));
	analyticframe->GetYaxis()->SetTitle("Counts");
	model_analytic.fitTo(hh, PrintLevel(-1));
	hh.plotOn(analyticframe);
	model_analytic.plotOn(analyticframe, LineColor(28));
	TString chiquadrato = TString::Format("#chi^{2} = %.1f, ndof = %d",analyticframe->chiSquare()*(30), (30-3));
	model_analytic.paramOn(analyticframe, Label(chiquadrato));
	//Visualize data
	//c[i] = new TCanvas(Form("d%d", i),Form("d%d", i), 800, 800);
	//analyticframe->Draw();
	
	std::cout << "Frequency " << i << " Number of entries: " << Hist->GetEntries() << std::endl;
	}

auto canvas1 = new TCanvas("1","1", 800,800);
auto g = new TGraph(vi.size(), vi.data(), Counts.data());
auto gres = new TGraphErrors(vi.size(), frequency.data(), resGas.data(), errorx.data(), errGas.data());
TString gtitle = TString::Format("Counts versus Frequencies");
g->SetTitle(gtitle);
g->SetMarkerColor(2);
g->SetMarkerStyle(23);
g->SetMinimum(-10);
g->GetYaxis()->SetTitle("Counts");
g->GetXaxis()->SetTitle("frequencies");
g->Draw("ap");
gres->SetMarkerStyle(23);
gres->SetMarkerColor(3);
gres->Draw("p");
}


void AllDataFit(){

//Generate dictionary
gInterpreter->GenerateDictionary("fitAllFrequency", "Headers/fitAllFrequencyh");
TString runNumber = "r68498"; // runlist r68465 r68481 r68489 r68498
std::vector<std::string> FileList;
FileList = getFiles(runNumber);

if (FileList.size() == 0) {
gROOT->ProcessLine(".x Conversion.cpp");
FileList = getFiles(runNumber);
}

double ExpectedRateMuons = (10.2)*FileList.size();

RooRealVar Ncosmic("Nfit_{cosmic}", "Nfit_{cosmic}" , ExpectedRateMuons, ExpectedRateMuons);
ROOT::RDataFrame rdf("myTree", FileList);
auto displ = rdf.Display({"CutsType0","CutsType1","CutsType2", "X", "Y", "Z"}, 5);
displ->Print();
auto rdf2 = rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto Hist = rdf2.Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius");


//Define variable to store the data
RooRealVar x("x", "r [cm]", 0, 4);
x.setBins(30);
//retrieve the histogram
RooDataHist hh("dh3", "dh3", x, Import(*Hist));
//RAYLEIGH
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay)))", RooArgSet(x,sigRay));
//COSMIC
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));
//MIXING
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("#sigma_{Mix}", "sigMix",0.85, 0.85,0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

//COMPLETE MODEL
RooRealVar Nmix("Nfit_{mix}","Nfit_{mix}",0.,2000.);
RooRealVar Ngas("Nfit_{gas}","Nfit_{gas}", -100,2000);
RooAddPdf model_analytic("model", "model", RooArgList{gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix,Ngas,Ncosmic});
//Do the fit
RooPlot *analyticframe = x.frame(Title("All the Data"));
analyticframe->GetYaxis()->SetTitle("Counts");
model_analytic.fitTo(hh, PrintLevel(-1));
hh.plotOn(analyticframe);
model_analytic.plotOn(analyticframe, LineColor(28));
TString chiquadrato = TString::Format("#chi^{2} = %.1f, ndof = %d",analyticframe->chiSquare()*(30), (30-3));
model_analytic.paramOn(analyticframe, Label(chiquadrato));
auto canvas = new TCanvas("Fit all the data together","d",800,800);
analyticframe->Draw();
}
