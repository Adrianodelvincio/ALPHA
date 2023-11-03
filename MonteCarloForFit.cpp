#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
using namespace RooFit;

void MonteCarloForFit(){
ROOT::RDataFrame mix("myTree", "Spectroscopy/RootCut1Data/MixCut1.root");
ROOT::RDataFrame uwlosses("myTree", "Spectroscopy/RootCut1Data/UWCut1.root");
ROOT::RDataFrame cosmic("myTree", "Spectroscopy/RootCut1Data/CosmicCut1.root");

RooRealVar x("x","x",0.,4.);

//RETRIEVE PDF MIXING
auto histMix = mix.Histo1D({"Mixing","Counts",30u,0.,4.}, "Radius");
RooDataHist mix_h("dh0", "dh0", x, Import(*histMix));
RooHistPdf PdfMixing("mixingpdf", "mixingpdf", x, mix_h, 0);

//PDF UWlosses
//FIXING sigRay to the value from the fit to the model (see analysisMLE.cpp)
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " x/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

//PDF Cosmic
//RooRealVar q("q", "q",0.25, 0.25 - 200);
RooGenericPdf linearFit("linearFit", "linear model", "(0.125)*x", RooArgSet(x));

//Create a PDF with as the sum of all three PDF

RooRealVar Nmix_t("Nmix","Nmix",100);
RooRealVar Nuw_a("Nuw","Nuw", 50);
RooRealVar Nbk_a("Nbk", "Nbk", 20);
RooPlot *frame1 = x.frame(Title("Sum Of Three Pdfs"));
RooAddPdf model_analytic("model", "model", RooArgList{PdfMixing,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});

RooPlot *frame2 = x.frame(Title("Sum of PdfMixing + rayleigh"));
RooAddPdf mix_cosmic("model1", "model1", RooArgList{PdfMixing, Rayleigh}, RooArgList{Nmix_t,Nuw_a});

RooPlot *frame3 = x.frame(Title("Sum of Rayleigh + cosmic"));
RooAddPdf rayleigh_cosmic("model2", "model2", RooArgList{Rayleigh,linearFit}, RooArgList{Nuw_a,Nbk_a});

RooPlot *frame4 = x.frame(Title("Mixing from the hist"));
RooAddPdf mix_basta("model3", "model3", RooArgList{PdfMixing}, RooArgList{Nmix_t});

model_analytic.plotOn(frame1, LineColor(28));
model_analytic.paramOn(frame1);
mix_cosmic.plotOn(frame2, LineColor(27));
mix_basta.plotOn(frame4,LineColor(26));
rayleigh_cosmic.plotOn(frame3, LineColor(26));

auto canvas1 = new TCanvas("d1","d1",800,800);
frame1->Draw();

auto canvas2 = new TCanvas("d2", "d2",800,800);
frame2->Draw();

auto canvas3 = new TCanvas("d3", "d3", 800,800);
frame3->Draw();

auto canvas4 = new TCanvas("d4","d4",800,800);
frame4->Draw();
}
