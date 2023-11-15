#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
#include "TGraphErrors.h"
using namespace RooFit;


void fitAllFrequency(){
// DATA TO FIT
TString cartella = TString::Format("Dataset/");
//TString firstpart = TString::Format("r68465_uw_exp_freq");
//TString firstpart = TString::Format("r68481_uw_exp_freq");
TString firstpart = TString::Format("r68481_uw_exp_freq");

TString formato = TString::Format(".csv");
vector<int> vi{1,2,3,4,5,6,7,8};
vector<double> frequency;
vector<int> Counts;
vector<double> resGas;
vector<double> errorx;
vector<double> errGas;
TCanvas *c[8];

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

for (int i : vi){
	TString secondpart = TString::Format("%d.vertex", i);
	frequency.push_back(i);
	//ROOTDATAFRAME MACHINERY
	auto rdf = ROOT::RDF::MakeCsvDataFrame(cartella + firstpart+ secondpart + formato);
	//auto displ = rdf.Display({"CutsType0","CutsType1","CutsType2", "X", "Y", "Z"}, 5);
    	//displ->Print();
	auto rdf2 = rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
	auto Hist = rdf2.Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius");
	//Save the entries
	Counts.push_back(Hist->GetEntries());
	resGas.push_back(static_cast<int> (Ngas.getVal()));
	errGas.push_back(static_cast<int> (Ngas.errorVar()->getVal()));
	errorx.push_back(0.);
	RooDataHist hh("dh3", "dh3", x, Import(*Hist));
	RooPlot *analyticframe = x.frame(Title(firstpart +  secondpart));
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

//ROOT::RDataFrame total_rdf("myTree", "DataSetROOT/r68465_cut1.root");
//ROOT::RDataFrame total_rdf("myTree", "DataSetROOT/r68481_cut1.root");
//ROOT::RDataFrame total_rdf("myTree", "DataSetROOT/r68489_cut1.root");
ROOT::RDataFrame total_rdf("myTree", "DataSetROOT/r68498_cut1.root");
RooRealVar Ncosmic("Nfit_{cosmic}", "Nfit_{cosmic}", 70, 70);

//Display some data
auto displ = total_rdf.Display({"CutsType0","CutsType1","CutsType2", "X", "Y", "Z"}, 5);
displ->Print();
auto rdf2 = total_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
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

