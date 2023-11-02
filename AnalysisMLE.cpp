#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
using namespace RooFit;

void AnalysisMLE(){

//MIXING
ROOT::RDataFrame mix_rdf("myTree", {"DataSetROOT/r68814_mixing.vertex.root",
		"DataSetROOT/r68839_mixing.vertex.root",
		"DataSetROOT/r68859_mixing.vertex.root",
		"DataSetROOT/r68871_mixing.vertex.root",
		"DataSetROOT/r68903_mixing.vertex.root",
		"DataSetROOT/r68905_mixing.vertex.root",
		"DataSetROOT/r68927_mixing.vertex.root",
		"DataSetROOT/r69126_mixing.vertex.root",
		"DataSetROOT/r69142_mixing.vertex.root"});

auto mix2_rdf = mix_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto histMix = mix2_rdf.Histo1D({"Mixing","Counts",30u,0.,4.}, "Radius");

RooRealVar x("x", "x", 0, 4);
x.setBins(30);

RooDataHist mix_h("dh0", "dh0", x, Import(*histMix));
RooPlot *frame1 = x.frame(Title("Mixing PDF"));
RooHistPdf PdfMixing("mixingpdf", "mixingpdf", x, mix_h, 0);
mix_h.plotOn(frame1);
PdfMixing.plotOn(frame1);

auto canvas0 = new TCanvas("d0", "d0", 800,800);
frame1->Draw();

// UWLOSSES
ROOT::RDataFrame uw_rdf("myTree",{"DataSetROOT/r68814_uwlosses_160.vertex.root",
	"DataSetROOT/r68839_uwlosses_160.vertex.root",
	"DataSetROOT/r68859_uwlosses_160.vertex.root",
	"DataSetROOT/r68871_uwlosses_160.vertex.root",
	"DataSetROOT/r68903_uwlosses_160.vertex.root",
	"DataSetROOT/r68905_uwlosses_160.vertex.root",
	"DataSetROOT/r68927_uwlosses_160.vertex.root",
	"DataSetROOT/r69126_uwlosses_160.vertex.root",
	"DataSetROOT/r69142_uwlosses_160.vertex.root"});

auto uw2_rdf = uw_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto histUw = uw2_rdf.Histo1D({"UWlosses","Counts",30u,0.,4.}, "Radius");

RooDataHist Uw_h("dh1", "dh1", x, Import(*histUw));
RooPlot *frame2 = x.frame(Title("UW losses PDF"));
RooHistPdf PdfUwlosses("uwlossespdf", "uwlossespdf", x, Uw_h, 0);
Uw_h.plotOn(frame2);
PdfUwlosses.plotOn(frame2);

auto canvas1 = new TCanvas("d1", "d1",800,800);
frame2->Draw();

//COSMIC

ROOT::RDataFrame cosmic_rdf("myTree",{"DataSetROOT/r68949_cosmics.vertex.root",
	"DataSetROOT/r69177_cosmics.vertex.root",
	"DataSetROOT/r69207_cosmics.vertex.root",
	"DataSetROOT/r69219_cosmics.vertexroot"});

auto cosmic2_rdf = cosmic_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto histCosmic = cosmic2_rdf.Histo1D({"Cosmic","Counts",30u,0.,4.}, "Radius");

RooDataHist cosm_h("dh2", "dh2", x, Import(*histCosmic));
RooPlot *frame3 = x.frame(Title("Cosmic PDF"));
RooHistPdf PdfBk("Cosmic", "cosmic", x, cosm_h, 0);
cosm_h.plotOn(frame3);
PdfBk.plotOn(frame3);

auto canvas3 = new TCanvas("d2","d2",800,800);
frame3->Draw();

// DATA TO FIT
auto f4_rdf = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Dataset/r68465_uw_exp_freq4.vertex.csv");
auto displ = f4_rdf.Display({"CutsType0","CutsType1","CutsType2", "X", "Y", "Z"}, 5);
    displ->Print();

auto f4_2_rdf = f4_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto Histf4_h = f4_2_rdf.Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius"); 


RooDataHist f4_h("dh3", "dh3", x, Import(*Histf4_h));


// Fit a Gaussian pdf to the data
RooRealVar mean("mean", "mean", 0, -10, 10);
RooRealVar sigma("sigma", "sigma", 3, 0.1, 10);
RooGaussian gauss("gauss", "gauss", x, mean, sigma);
gauss.fitTo(f4_h, PrintLevel(-1));
RooPlot *frameG = x.frame(Title("fit with gaussian"));
f4_h.plotOn(frameG);
gauss.plotOn(frameG);
gauss.paramOn(frameG);

auto canvas5 = new TCanvas("d5", "d5",800,800);
frameG->Draw();

RooRealVar Nmix("Nmix", "Nmix", 0., 200);
RooRealVar Nuw("Nuw", "Nuw", 0., 200);
RooRealVar Nbk("Nbk", "Nbk", 0., 200);
RooPlot *frame4 = x.frame(Title("Fit to f4 Dataset"));
RooAddPdf model("model","model", RooArgList{PdfMixing,PdfUwlosses,PdfBk}, RooArgList{Nmix, Nuw, Nbk} );
model.fitTo(f4_h,PrintLevel(-1));
f4_h.plotOn(frame4);
model.plotOn(frame4);
auto canvas4 = new TCanvas("d4","d4", 800, 800);
frame4->Draw();
}
