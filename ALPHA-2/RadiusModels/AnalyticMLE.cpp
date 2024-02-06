#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
using namespace RooFit;

void AnalyticMLE( TString filename = "r68465_f5.root"
				){

// DATA TO FIT
TString cartella = TString::Format("DataSetROOT/");

RooRealVar x("x", "r [cm]", 0, 4);
x.setBins(30);

// DISTRIBUTIONS FOR THE MLE FIT
// Residual Gas
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay)))", RooArgSet(x,sigRay));
// Cosmic
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));
// Mixing
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

// MODEL = Gaussian(mix) + Rayleigh(res. gas) + linear(cosmic)
RooRealVar Nmix_t("Nfit_{mix}","Nfit_{mix}",0.,200.);
RooRealVar Nuw_a("Nfit_{gas}","Nfit_{gas}", -100,100);
RooRealVar Nbk_a("Nfit_{cosmic}", "Nfit_{cosmic}",10.2,10.2);
RooAddPdf model_analytic("model", "model", RooArgList{gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});

// LOAD the data
ROOT::RDataFrame f4_rdf("myTree" ,{cartella + filename});
auto displ = f4_rdf.Display({"CutsType0","CutsType1","CutsType2", "X", "Y", "Z"}, 5);
displ->Print(); // Print on the shell some rows

auto f4_2_rdf = f4_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\""); //Cut 1 selection
auto Histf = f4_2_rdf.Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius"); 
RooDataHist hf("dh3", "dh3", x, Import(*Histf)); 


RooPlot *analyticframe = x.frame(Title(filename));
analyticframe->GetYaxis()->SetTitle("Counts");
model_analytic.fitTo(hf, PrintLevel(-1));
hf.plotOn(analyticframe);
model_analytic.plotOn(analyticframe, LineColor(28));
TString chiquadrato = TString::Format("#chi^{2} = %.1f, ndof = %d",analyticframe->chiSquare()*(30), (30-3));
model_analytic.paramOn(analyticframe, Label(chiquadrato),  Layout(0.12,0.5,0.9));

auto canvas = new TCanvas("MLE FIT","MLE FIT",800,800);
analyticframe->getAttText()->SetTextSize(0.03);
analyticframe->Draw();
}
