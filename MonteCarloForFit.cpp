#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
#include "TGraphErrors.h"
using namespace RooFit;


void SetVectors(int i,std::vector<double> &,std::vector<double> &,std::vector<double> &, RooRealVar,RooRealVar,RooRealVar,RooRealVar,RooRealVar, RooRealVar);
Double_t media(std::vector<Double_t>&v);


void MonteCarloForFit(){
double a = 0.33, b = 0.33, c = 0.33; // percentage of the Pdfs
// N: statistic of the data, Nloop number of iterations
int N = 1000, Nloop = 1000; // f4 contains 165 event, use this number for real simulation
//if the number of cosmic is fixed (costant time window), then use Nfix for Nbk_a
//For the cosmic realistic scenario
double Nfix = 0;
//Nfix = 10.2;
//c = Nfix/N
//a = a - c/2;
//b = b - c/2;

//CREATE AN ISTOGRAM TO SAVE THE QUANTITIES OF THE MONTECARLO
TH1* hmix = new TH1I("h1", "MIX : (N - Nfit)/sigma",30, -6,6);
TH1* huw =  new TH1I("h2", "UW : (N - Nfit)/sigma",30,-6,6);
TH1* hbk =  new TH1I("h3", "BK : (N - Nfit)/sigma",30,-6,6);
TH1* hchisq = new TH1I("hchisq", "chiquadro", 30, 0, 50);

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

//MIXING analytic model
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

RooRealVar Nmix_t("Nmix","Nmix",a*N, -3000, +3000);
RooRealVar Nuw_a ("Nuw", "Nuw", b*N, -3000, +3000);
RooRealVar Nbk_a ("Nbk", "Nbk", c*N, -3000, +3000);
RooPlot *frame1 = x.frame(Title("Sum Of Three Pdfs"));
//Model to generate the data
RooAddPdf model_gen("model", "model", RooArgList{/*PdfMixing*/ gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});
model_gen.plotOn(frame1, LineColor(28));
model_gen.paramOn(frame1);

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
auto canvass1 = new TCanvas("d2", "d2",800,800);
frame2->Draw();

auto canvass2 = new TCanvas("d3", "d3", 800,800);
frame3->Draw();

auto canvass3 = new TCanvas("d4","d4",800,800);
frame4->Draw();
*/

//Generate a dataset of N events in x from model_analytic pdf
std::unique_ptr<RooDataSet> data{model_gen.generate(x,N)}; //generate the data
RooDataHist *histXgen = data->binnedClone(); //create a binned histogram
RooPlot *frame2 = x.frame(Title("Toy Model"));
//data->plotOn(frame2);
//model_gen.plotOn(frame2);
//model_gen.paramOn(frame2);
histXgen->plotOn(frame2);

// FIT the toy model
//variable of the fit
RooRealVar Nmix_f("Nmixf","Nmix",a*N,-3000, +3000);
RooRealVar Nuw_f("Nuwf","Nuw",b*N,-3000, +3000);
RooRealVar Nbk_f("Nbkf", "Nbk",c*N,-3000,3000);
//create the fit model
RooAddPdf model_forfit("model2", "model", RooArgList{/*PdfMixing*/ gauss_Mix,Rayleigh,linearFit},RooArgList{Nmix_f,Nuw_f,Nbk_f});
RooPlot *frame3 = x.frame(Title("Fit Toy Model"));
std::unique_ptr<RooFitResult> fitResult(model_forfit.fitTo(*data, Save()));
histXgen->plotOn(frame3);
model_forfit.plotOn(frame3,LineColor(27));

//TString cosmici = TString::Format("Cosmici Nbk = %.2f ", Nfix*N);
model_forfit.paramOn(frame3/*,Label(cosmici)*/);
double fit0, fit1, fit2;
fit0 = Nmix_f.getVal();
fit1 = Nuw_f.getVal();
fit2 = Nbk_f.getVal();

//Print expected number of event from the fit
//fitResult->Print();
std::cout << "\n \n" << std::endl;
std::cout << "Nmix, Nuw, Nbk expected from the tow" << std::endl;
std::cout << "Nmix: " << N*a << " Nuw: " << N*b << " Nbk: " << c*N << std::endl;
std::cout << "values from the fit " << std::endl;
std::cout << "Nmix: " << fit0 << " Nuw: " << fit1 << " Nbk " << fit2 << std::endl;

double fit3 = 0,fit4 = 0,fit5 = 0;
//Disable some unuseful print
RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
//LOOP ON THE FIT NOW
for(int i = 0; i < Nloop; i++){
	//generate the data
	std::unique_ptr<RooDataSet> dataLoop{model_gen.generate(x,N)};
	model_forfit.fitTo(*dataLoop, Save(), PrintLevel(-1));
	//save and store the fit values
	fit3 = (Nmix_f.getVal() - N*a)/Nmix_f.errorVar()->getVal();
	fit4 = (Nuw_f.getVal() - N*b)/Nuw_f.errorVar()->getVal();
	fit5 = (Nbk_f.getVal() - N*c)/Nbk_f.errorVar()->getVal();
	//Get the ChiSquare
	RooPlot *frameLoop = x.frame(Title("data"));
	dataLoop->plotOn(frameLoop);
	model_forfit.plotOn(frameLoop);
	//fill histogram
	hmix->Fill(fit3);
	huw->Fill(fit4);
	hbk->Fill(fit5);
	hchisq->Fill((frameLoop->chiSquare())*(30));
	if(i%10 == 0){
		std::cout << "Event LOOP N°" << i << std::endl;
		std::cout << "fit:      " << "Nmix: " << Nmix_f.getVal() << " Nuw: " << Nuw_f.getVal() << " Nbk " << Nbk_f.getVal() << std::endl;
		std::cout << "Expected: " << "Nmix: " << N*a << " Nuw: " << N*b << " Nbk: " << Nfix*N << std::endl;
	}
}

//Fit the Histograms of N - Ntrue
hmix->Fit("gaus");
huw->Fit("gaus");
hbk->Fit("gaus");


auto legend = new TLegend(0.1,0.7,0.48,0.9);
legend->SetHeader("PDF Composition","C"); // option "C" allows to center the header
TString coeffMix = TString::Format("degree mixing = %.2f %% ", a*100);
TString coeffUw = TString::Format("degree uw = %.2f %%", b*100);
TString coeffbk = TString::Format("degree bk = %.2f %%", c*100);
legend->AddEntry(hmix,coeffMix);
legend->AddEntry(hmix,coeffUw);
legend->AddEntry(hmix,coeffbk);
auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
legend1->SetHeader("PDF Composition","C"); // option "C" allows to center the header
legend1->AddEntry(huw,coeffMix);
legend1->AddEntry(huw,coeffUw);
legend1->AddEntry(huw,coeffbk);
auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
legend2->SetHeader("PDF Composition","C"); // option "C" allows to center the header
legend2->AddEntry(hbk,coeffMix);
legend2->AddEntry(hbk,coeffUw);
legend2->AddEntry(hbk,coeffbk);


auto canvas1 = new TCanvas("d1","d1",800,800);
frame1->Draw();

auto canvas2 = new TCanvas("d2", "Toy Model data",800,800);
frame2->Draw();

auto canvas3 = new TCanvas("d3", "Fit to toy model", 800, 800);
frame3->Draw();

auto canvas4 = new TCanvas("d4", "Toy Result", 1000,1000);
auto pad = new TPad("pad", "pad",0,0,1,1);
gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
pad->Divide(3,1,0.,0.);
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
//legend->Draw();
pad->cd(3);
hbk->SetLineColor(1);
hbk->SetLineWidth(2);
hbk->Draw();

auto canvas5 = new TCanvas("d5","d5",800,800);
hchisq->Draw();
}



void SecondMontecarlo(){
double a = 0.90, b = 0.10, c = 0.10; // percentage of the Pdfs
// N: statistic of the data, Npoint number of iterations
int N = 500, Npoint = 30; // f4 contains 165 event, use this number for real simulation
int Nnested = 100;
double Nfix = 0;
//Nfix = 10.2;
//c = Nfix/N
//a = a - c/2;
//b = b - c/2;

RooRealVar x("x","x",0.,4.);
x.setBins(30);

//MIXING analytic model
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

//Model to generate the data
RooRealVar Nmix_t("Nmix","Nmix",a*N, -(2*N), +(2*N));
RooRealVar Nuw_a ("Nuw", "Nuw", b*N, -(2*N), +(2*N));
RooRealVar Nbk_a ("Nbk", "Nbk", c*N, -(2*N), +(2*N));
RooAddPdf model_gen("model", "model", RooArgList{/*PdfMixing*/ gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});

//create the fit model
RooRealVar Nmix_f("Nmixf","Nmix",a*N,-(2*N), +(2*N));
RooRealVar Nuw_f ("Nuwf", "Nuw", b*N,-(2*N), +(2*N));
RooRealVar Nbk_f ("Nbkf", "Nbk", c*N,-(2*N), +(2*N));
RooAddPdf model_forfit("model2", "model", RooArgList{/*PdfMixing*/ gauss_Mix,Rayleigh,linearFit},RooArgList{Nmix_f,Nuw_f,Nbk_f});

//Arrays to store the data
vector<double> mix;
vector<double> uw;
vector<double> bk;
vector<double> mixErrors;
vector<double> weight;
vector<double> Errors;

vector<double> mix_l;
vector<double> uw_l;
vector<double> bk_l;
vector<double> mixErrors_l;
vector<double> weight_l;

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
for(int i = 0; i < Npoint; i++){
	//Fix the weight of the pdf
	double wmix = c + static_cast<double>(i)*(1 - 2*c)/Npoint;
	weight.push_back(wmix);
	//Fix Pdf weights
	Nmix_t.setVal(N*wmix);
	Nuw_a.setVal(( 1-wmix-c)*N); //Expected event for uw is 1 - wmix - wbk
	Nbk_a.setVal(N*c);
	
	
	for(int j = 0; j < Nnested; j++){
		//generate the data
		std::unique_ptr<RooDataSet> dataLoop{model_gen.generate(x,N)};
		model_forfit.fitTo(*dataLoop, Save(), PrintLevel(-1));
		//save the data
		//set parameter of the fit
		SetVectors(j,mix_l,uw_l,bk_l,Nmix_f,Nuw_f,Nbk_f,Nmix_t,Nuw_a,Nbk_a);
		//set error of the parameters
		mixErrors_l.push_back(Nmix_f.errorVar()->getVal());
	}
	//save and store the fit values
	mix.push_back(media(mix_l));
	uw.push_back(media(uw_l));
	bk.push_back(media(bk_l));
	mixErrors.push_back(media(mixErrors_l));
	Errors.push_back(media(mixErrors_l)/sqrt(Nnested));
	//std::cout << media(mixErrors_l)/Nnested << std::endl;
	mix_l.clear();
	uw_l.clear();
	bk_l.clear();
	mixErrors_l.clear();
	if(i%1 == 0){
	std::cout << "Event LOOP N°" << i << std::endl;
	std::cout << "weight Mix: " << wmix << std::endl;
	std::cout << "fit:      " << "Nmix: " << mix[i] << " Nuw: " << uw[i] << " Nbk " << bk[i] << std::endl;
	std::cout << "Expected: " << "Nmix: " << Nmix_t.getVal() << " Nuw: " << Nuw_a.getVal() << " Nbk: " << Nbk_a.getVal() << std::endl;
	}
}

//Visualize the Montecarlo
vector<double> errorx(30,0);
auto g = new TGraphErrors(weight.size(), weight.data(), mix.data(),errorx.data(),Errors.data());

g->SetTitle("#mu_{N - Nfit}");
g->GetYaxis()->SetTitle("#mu_{N - Nfit}");
g->GetXaxis()->SetTitle("weight");
g->SetMinimum(-30);
g->SetMaximum(+30);
g->GetYaxis()->SetLimits(-20,20);
g->SetMarkerStyle(21);

auto canvas0 = new TCanvas("canvas0", "Nmix vs weight",1000,1000);
g->Draw("ap");

auto g1 = new TGraph(weight.size(),weight.data(),mixErrors.data());
g1->SetMarkerStyle(21);
g1->GetYaxis()->SetTitle("#sigma_{Nmix} vs. weight");
g1->SetMinimum(0.);
g1->SetMaximum(100.);
g1->GetYaxis()->SetLimits(-50,50);
//g1->SetMarkerStyle(21);

auto canvas1 = new TCanvas("canvas1", "#sigma_{mix} vs weight",1000,1000);
g1->Draw("ap");
}




void SetVectors(int i,std::vector<double> &mix,std::vector<double> &uw, std::vector<double> &bk, RooRealVar a,RooRealVar b,RooRealVar c,RooRealVar aa,RooRealVar bb,RooRealVar cc){
		mix.push_back(a.getVal() - aa.getVal());
		uw.push_back(b.getVal() - bb.getVal());
		bk.push_back(c.getVal() - cc.getVal());}
		

Double_t media(std::vector<Double_t>&v) {
if(v.size()<=0.) return 0.;
return (Double_t)std::accumulate(v.begin(), v.end(), 0.0) / (Double_t)v.size();
}
