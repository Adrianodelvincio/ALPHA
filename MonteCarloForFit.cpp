#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
using namespace RooFit;

void MonteCarloForFit(){

double a = 50, b = 50; // percentage of the Pdfs
int N = 165; // Number of event to generate
int Nloop = 2000; //loop lenght
//if the number of cosmic is fixed (costant time window), then use Nfix for Nbk_a
double Nfix = 0;
//Nfix = 10.2/N;

//CREATE AN ISTOGRAM TO SAVE THE QUANTITIES OF THE MONTECARLO
TH1* hmix = new TH1I("h1", "Nmix - Nmix_fit",30,-N-200,N+200);
TH1* huw = new TH1I("h2", "Nuw - Nuw_fit",30,-N-200,N+200);

//LOAD THE DATA
ROOT::RDataFrame mix("myTree", "Spectroscopy/RootCut1Data/MixCut1.root");
ROOT::RDataFrame uwlosses("myTree", "Spectroscopy/RootCut1Data/UWCut1.root");
ROOT::RDataFrame cosmic("myTree", "Spectroscopy/RootCut1Data/CosmicCut1.root");

RooRealVar x("x","x",0.,4.);
x.setBins(30);

//RETRIEVE PDF MIXING
auto histMix = mix.Histo1D({"Mixing","Counts",30u,0.,4.}, "Radius");
RooDataHist mix_h("dh0", "dh0", x, Import(*histMix));
RooHistPdf PdfMixing("mixingpdf", "mixingpdf", x, mix_h, 0);

//for mixing we can use also a gaussian distribution, that works fine
RooRealVar mu("mu", "mu", 2.38,2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85,0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);


//PDF UWlosses
//FIXING sigRay to the value from the fit to the model (see analysisMLE.cpp)
RooRealVar sigRay("sigRay", "sigma", 1.722,1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

//PDF Cosmic
//RooRealVar q("q", "q",0.25, 0.25 - 200);
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));

//Create a PDF with as the sum of all three PDF
// Pay attention, the first value of Nmix,Nuw,Nbk are used to generate the data, and represent the percentage of the pdf that are used inside the final model.
a = a*(N - Nfix*N)/N;
b = b*(N - Nfix*N)/N;

RooRealVar Nmix_t("Nmix","Nmix", a, -3000, +3000);
RooRealVar Nuw_a("Nuw","Nuw", b, -3000, +3000);
RooRealVar Nbk_a("Nbk", "Nbk", Nfix);
RooPlot *frame1 = x.frame(Title("Sum Of Three Pdfs"));
RooAddPdf model_analytic("model", "model", RooArgList{PdfMixing /*gauss_Mix*/,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});
model_analytic.plotOn(frame1, LineColor(28));
model_analytic.paramOn(frame1);

/* SINGLE DISTRIBUTIONS
RooPlot *frame2 = x.frame(Title("Sum of PdfMixing + rayleigh"));
RooAddPdf mix_cosmic("model1", "model1", RooArgList{PdfMixing, linearFit}, RooArgList{Nmix_t,Nuw_a});

RooPlot *frame3 = x.frame(Title("Sum of Rayleigh + cosmic"));
RooAddPdf rayleigh_cosmic("model2", "model2", RooArgList{Rayleigh,linearFit}, RooArgList{Nuw_a,Nbk_a});

RooPlot *frame4 = x.frame(Title("Mixing from the hist"));
RooAddPdf mix_basta("model3", "model3", RooArgList{PdfMixing}, RooArgList{Nmix_t});

mix_cosmic.plotOn(frame2, LineColor(27));
mix_basta.plotOn(frame4,LineColor(26));
rayleigh_cosmic.plotOn(frame3, LineColor(26));

// SINGLE DISTRIBUTIONS
auto canvas2 = new TCanvas("d2", "d2",800,800);
frame2->Draw();

auto canvas3 = new TCanvas("d3", "d3", 800,800);
frame3->Draw();

auto canvas4 = new TCanvas("d4","d4",800,800);
frame4->Draw();
*/




//Generate a dataset of N events in x from model_analytic pdf
std::unique_ptr<RooDataSet> data{model_analytic.generate(x,N)}; //generate the data
RooDataHist *histXgen = data->binnedClone(); //create a binned histogram
RooPlot *frame2 = x.frame(Title("Toy Model"));
histXgen->plotOn(frame2);
//data->plotOn(frame2);
model_analytic.plotOn(frame2);
model_analytic.paramOn(frame2);

// FIT the toy model
RooAddPdf model_forfit("model2", "model", RooArgList{PdfMixing /*gauss_Mix*/ ,Rayleigh,linearFit},RooArgList{Nmix_t,Nuw_a,Nbk_a});
RooPlot *frame3 = x.frame(Title("Fit Toy Model"));
std::unique_ptr<RooFitResult> fitResult(model_forfit.fitTo(*data, Save()));
histXgen->plotOn(frame3);
model_forfit.plotOn(frame3,LineColor(27));
TString cosmici = TString::Format("Cosmici Nbk = %.2f ", Nfix*N);
model_forfit.paramOn(frame3, Label(cosmici));
double fit0, fit1, fit2;
fit0 = Nmix_t.getVal();
fit1 = Nuw_a.getVal();
fit2 = Nbk_a.getVal();

double fit3 = 0,fit4 = 0,fit5 = 0;

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
//LOOP ON THE FIT NOW
for(int i = 0; i < Nloop; i++){
	std::cout << "Event LOOP N°" << i << std::endl; 
	//generate the data
	std::unique_ptr<RooDataSet> dataLoop{model_analytic.generate(x,N)};
	model_forfit.fitTo(*dataLoop, Save(), PrintLevel(-1));
	//save and store the fit values
	fit3 = (Nmix_t.getVal() - N*a/100) /*/Nmix_t.errorVar()->getVal()*/;
	fit4 = (Nuw_a.getVal() - N*b/100)  /*/Nuw_a.errorVar()->getVal() */;
	fit5 = Nbk_a.getVal() -Nfix*N;
	//fill histogram
	hmix->Fill(fit3);
	huw->Fill(fit4);
	
	std::cout << "fit:      " << "Nmix: " << Nmix_t.getVal() << " Nuw: " << Nuw_a.getVal() << " Nbk " << Nbk_a.getVal() << std::endl;
	std::cout << "Expected: " << "Nmix: " << N*a/100 << " Nuw: " << N*b/100 << " Nbk: " << Nfix*N << std::endl;
}



//Print expected number of event from the fit
//fitResult->Print();
std::cout << "\n \n" << std::endl;
std::cout << "Nmix, Nuw, Nbk expected from the tow" << std::endl;
std::cout << "Nmix: " << N*a/100 << " Nuw: " << N*b/100 << " Nbk: " << Nfix*N << std::endl;
std::cout << "values from the fit " << std::endl;
std::cout << "Nmix: " << fit0 << " Nuw: " << fit1 << " Nbk " << fit2 << std::endl;

auto legend = new TLegend(0.1,0.7,0.48,0.9);
legend->SetHeader("PDF Composition","C"); // option "C" allows to center the header
TString coeffMix = TString::Format("degree mixing = %.2f %% ", a);
TString coeffUw = TString::Format("degree uw = %.2f %%", b);
TString coeffbk = TString::Format("degree bk = %.2f %%", Nfix);
legend->AddEntry(hmix,coeffMix);
legend->AddEntry(hmix,coeffUw);
legend->AddEntry(hmix,coeffbk);




auto canvas1 = new TCanvas("d1","d1",800,800);
frame1->Draw();

auto canvas2 = new TCanvas("d2", "Toy Model data",800,800);
frame2->Draw();

auto canvas3 = new TCanvas("d3", "Fit to toy model", 800, 800);
frame3->Draw();

auto canvas4 = new TCanvas("d4", "Toy Result", 800,800);
auto pad = new TPad("pad", "pad",0,0,1,1);
pad->Divide(2,1,0.02,0.02);
pad->Draw();

pad->cd(1);
hmix->SetLineColor(1);
hmix->SetLineWidth(2);
hmix->Draw();
legend->Draw();
pad->cd(2);
huw->SetLineColor(1);
huw->SetLineWidth(2);
huw->Draw();
legend->Draw();
}
