#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
#include "TGraphErrors.h"
#include "Headers/MontecarloForFit.h"

using namespace RooFit;

void MonteCarloForFit(){
// percentage of the Pdfs
double a = 0.33;
double b = 0.33;
double c = 0.33; 
int CosmicFixed = 0;
int N = 1000;
int Nloop = 1000; 
double Nfix = 10.2;
if(CosmicFixed == 1){
SetProb(&a, &b, &c, N, Nfix);
}

//Generate dictionary
gInterpreter->GenerateDictionary("Functions", "Headers/MontecarloForFit.h");

//CREATE AN ISTOGRAM TO SAVE THE QUANTITIES OF THE MONTECARLO
TH1* hmix = new TH1I("h1", "MIX : (N_{fit} - N_{gen})/#sigma",30, -6,6);
TH1* huw =  new TH1I("h2", "Residual gas : (N_{fit} - N_{gen})/#sigma",30,-6,6);
TH1* hbk =  new TH1I("h3", "Cosmic : (N_{fit} - N_{gen})/#sigma",30,-6,6);
TH1* hchisq = new TH1I("hchisq", "chiquadro", 30, 0, 50);

RooRealVar x("x","r [cm]",0.,4.);
x.setBins(30);

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

RooRealVar Nmix_t("Nmix","Nmix",a*N, -3000, +3000);
RooRealVar Nuw_a ("Ngas", "Ngas", b*N, -3000, +3000);
RooRealVar Nbk_a ("Ncosmic", "Ncosmic", c*N, -3000, +3000);
RooPlot *frame1 = x.frame(Title("Sum Of Three Pdfs"));
//Model to generate the data
RooAddPdf model_gen("model", "model", RooArgList{gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});
model_gen.plotOn(frame1, LineColor(28));
model_gen.paramOn(frame1);

//Generate a dataset of N events in x from model_analytic pdf
std::unique_ptr<RooDataSet> data{model_gen.generate(x,Extended())}; //generate the data
RooDataHist *histXgen = data->binnedClone(); //create a binned histogram
RooPlot *frame2 = x.frame(Title("Toy Model"));
//data->plotOn(frame2);
//model_gen.plotOn(frame2);
//model_gen.paramOn(frame2);
histXgen->plotOn(frame2);

// FIT the toy model
//variable of the fit
RooRealVar Nmix_f("Nfit_{mix}","Nmix",a*N,-3000, +3000);
RooRealVar Nuw_f("Nfit_{gas}","Ngas",b*N,-3000, +3000);
//RooRealVar Nuw_f("Nuwf","Nuw", 0.);
RooRealVar Nbk_f("Nfit_{cosmic}", "Nbk",c*N,-3000,3000);
//create the fit model
RooAddPdf model_forfit("model2", "model", RooArgList{gauss_Mix,Rayleigh,linearFit},RooArgList{Nmix_f,Nuw_f,Nbk_f});
RooPlot *frame3 = x.frame(Title("Toy Model Fit"));
frame3->GetYaxis()->SetTitle("Counts");
model_forfit.fitTo(*data, Save());
histXgen->plotOn(frame3);
model_forfit.plotOn(frame3,LineColor(27));
TString chiquadrato = TString::Format("#chi^{2} = %.1f, ndof = %d",frame3->chiSquare()*(30), (30-3));
//TString cosmici = TString::Format("Cosmici Nbk = %.2f ", Nfix*N);
model_forfit.paramOn(frame3,Label(chiquadrato));
//frame3->getAttText()->SetTextSize(9); 
double fit0, fit1, fit2;
fit0 = Nmix_f.getVal();
fit1 = Nuw_f.getVal();
fit2 = Nbk_f.getVal();

//Print expected number of event from the fit
//fitResult->Print();
std::cout << "\n \n" << std::endl;
std::cout << "Sample N = " << data->numEntries() << std::endl;
std::cout << "Nmix, Nuw, Nbk expected from the tow" << std::endl;
std::cout << "Nmix: " << N*a << " Nuw: " << N*b << " Nbk: " << c*N << std::endl;
std::cout << "values from the fit " << std::endl;
std::cout << "Nmix: " << fit0 << " Nuw: " << fit1 << " Nbk " << fit2 << std::endl;

double fit3 = 0,fit4 = 0,fit5 = 0;
//Disable some unuseful print
RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
TRandom *sampleN = new TRandom();
//LOOP ON THE FIT NOW
for(int i = 0; i < Nloop; i++){
	//generate the data
	std::unique_ptr<RooDataSet> dataLoop{model_gen.generate(x,Extended())};
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
		std::cout << "Event LOOP N°" << i << " Gen Event: " << dataLoop->numEntries() <<std::endl;
		std::cout << "fit:      " << "Nmix: " << Nmix_f.getVal() << " Nuw: " << Nuw_f.getVal() << " Nbk " << Nbk_f.getVal() << std::endl;
		std::cout << "Expected: " << "Nmix: " << N*a << " Nuw: " << N*b << " Nbk: " << c*N << std::endl;
	}
}

//Fit the Histograms of N - Ntrue
hmix->Fit("gaus");
huw->Fit("gaus");
hbk->Fit("gaus");


auto legend = new TLegend(0.1,0.7,0.48,0.9);
legend->SetHeader("PDF Composition","C"); // option "C" allows to center the header
TString coeffMix = TString::Format("degree mixing = %.2f %% ", a*100);
TString coeffUw = TString::Format("degree gas = %.2f %%", b*100);
TString coeffbk = TString::Format("degree cosmic = %.2f %%", c*100);
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


//auto canvas1 = new TCanvas("d1","d1",800,800);
//frame1->Draw();

auto canvas2 = new TCanvas("d2", "Toy Model data",800,800);
frame2->Draw();

auto canvas3 = new TCanvas("d3", "Fit to toy model", 800, 800);
TString filenamefit = TString::Format("PlotMLEfit/N%d/FitToy(%d,%d,%d).pdf",N,static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
frame3->Draw();
canvas3->SaveAs(filenamefit);

auto canvas4 = new TCanvas("d4", "Toy Result", 1000,550);
auto pad = new TPad("pad", "pad",0,0,1,1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
pad->Divide(3,1,0.,0.);
pad->Draw();
pad->cd(1);
hmix->SetLineColor(1);
hmix->SetLineWidth(2);
hmix->Draw();
//legend->Draw();
pad->cd(2);
huw->SetLineColor(1);
huw->SetLineWidth(2);
huw->Draw();
//legend->Draw();
pad->cd(3);
hbk->SetLineColor(1);
hbk->SetLineWidth(2);
hbk->Draw();
TString percorso = TString::Format("PlotMLEfit/N%d/",N);

if(CosmicFixed == 1){
TString nameFile0 = TString::Format("ToyNmix_with_bkFixed(%d,%d,%d).pdf",static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
//canvas4->SaveAs(percorso + nameFile0);
}
else{
TString nameFile0 = TString::Format("ToyNmix(%d,%d,%d).pdf",static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
//canvas4->SaveAs(percorso + nameFile0);
}
auto canvas5 = new TCanvas("d5","d5",800,800);
hchisq->Draw();

auto canvas6 = new TCanvas("d6", "Toy Result", 1000,700);
auto pad1 = new TPad("pad", "pad",0,0,1,1);
gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
pad1->Divide(2,1,0.,0.);
pad1->Draw();
pad1->cd(1);
hmix->SetLineColor(1);
hmix->SetLineWidth(2);
hmix->Draw();
legend->Draw();
pad1->cd(2);
hbk->SetLineColor(1);
hbk->SetLineWidth(2);
hbk->Draw();

//legend->Draw();
if(CosmicFixed == 1){
TString nameFile1 = TString::Format("ToyNmixNbk_with_bkFixed(%d,%d,%d).pdf",static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
canvas6->SaveAs(percorso + nameFile1);
}
else{
TString nameFile1 = TString::Format("ToyNmixNbk(%d,%d,%d).pdf",static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
canvas6->SaveAs(percorso + nameFile1);
	}
}

void SecondMontecarlo(){
// percentage of the Pdfs
double a = 0.9;
double b = 0. ;
double c = 0.10; 
// N: statistic of the data, Npoint number of iterations
int CosmicFixed = 0;
int N = 1000; 
int Npoint = 30; // f4 contains 165 event, use this number for real simulation
int Nnested = 100;

//Generate dictionary
gInterpreter->GenerateDictionary("MontecarloForFit", "Headers/MontecarloForFit.h");

// Fix the cosmic Background
double Nfix = 10.2;
if(CosmicFixed == 1){
SetProb(&a,&b,&c,N,Nfix);
}

// Set the variable to store the data
RooRealVar x("x","r [cm]",0.,4.);
x.setBins(30);

RooRealVar mu("mu", "mu", 2.38,2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85,0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

RooRealVar sigRay("sigRay", "#sigma_{rayleigh}", 1.722,1.722); // fixing sigma from the fit value
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x)); //PDF Cosmic

//Model to generate the data
RooRealVar Nmix_t("Nmix",    "Nmix", 	0, -(2*N), +(2*N));
RooRealVar Nuw_a ("Ngas",    "Ngas", 	0, -(2*N), +(2*N));
RooRealVar Nbk_a ("Ncosmic", "Ncosmic", 0, -(2*N), +(2*N));
RooAddPdf model_gen("model", "model", RooArgList{ gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});

//Model for the Fit
RooRealVar Nmix_f("Nmix fit",    "Nmix",    0,-(2*N), +(2*N));
RooRealVar Nuw_f ("Ngas fit",    "Ngas",    0,-(2*N), +(2*N));
RooRealVar Nbk_f ("Ncosmic fit", "Ncosmic", 0,-(2*N), +(2*N));
RooAddPdf model_forfit("model2", "model", RooArgList{gauss_Mix,Rayleigh,linearFit},RooArgList{Nmix_f,Nuw_f,Nbk_f});

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
	double wmix = 0 + static_cast<double>(i)*(1 - c)/(Npoint-1);
	weight.push_back(wmix);
	ChangeWeight(&Nmix_t, &Nuw_a, &Nbk_a, wmix, c, N); // Change the weight
	for(int j = 0; j < Nnested; j++){
		//generate the data
		std::unique_ptr<RooDataSet> dataLoop{model_gen.generate(x,Extended())};
		if(j%10 == 0 ){
		std::cout << "Event gen " << dataLoop->numEntries() << std::endl;}
		model_forfit.fitTo(*dataLoop, Save(), PrintLevel(-1));
		SetVectors(j,mix_l,uw_l,bk_l,Nmix_f,Nuw_f,Nbk_f,Nmix_t,Nuw_a,Nbk_a); //set parameter of the fit
		mixErrors_l.push_back(Nmix_f.errorVar()->getVal()); //set error of the parameters
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
	PrintInfo(Nmix_f,Nuw_f,Nbk_f,Nmix_t, Nuw_a, Nbk_a, wmix, i);
	}
}

//Visualize the Montecarlo
vector<double> errorx(30,0);
auto g = new TGraphErrors(weight.size(), weight.data(), mix.data(),errorx.data(),Errors.data());
auto cosmic_g = new TGraph(weight.size(), weight.data(),bk.data());
std::ostringstream s; s << "MIX: N_{gen} - N_{fit} averaged over " << Nnested << " trials";
g->SetTitle(s.str().c_str());
g->GetYaxis()->SetTitle("<N_{gen} - N_{fit}>");
g->GetXaxis()->SetTitle("a");

if(N > 500){
	g->SetMinimum(-50);
	g->SetMaximum(50);
}
else{
	g->SetMinimum(-30);
	g->SetMaximum(+30);
}

g->GetYaxis()->SetLimits(-30,30);
g->SetMarkerStyle(21);
cosmic_g->SetMarkerStyle(20);
cosmic_g->SetMarkerColor(2);

auto canvas0 = new TCanvas("canvas0", "Nmix vs weight",1000,1000);
g->Draw("ap");
cosmic_g->Draw("p");

TString percorso = TString::Format("PlotMLEfit/N%d/",N);
TString nameFile0 = TString::Format("Nmix_minus_Nfit(%d,%d,%d).pdf",static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
canvas0->SaveAs(percorso + nameFile0);

auto g1 = new TGraph(weight.size(),weight.data(),mixErrors.data());
g1->SetMarkerStyle(21);
g1->SetTitle("MIX: #sigma vs weight");
g1->GetYaxis()->SetTitle("#sigma_{Nmix} vs. a");
g1->GetXaxis()->SetTitle("a");
if(N > 500){
g1->SetMinimum(50);
g1->SetMaximum(80);}
else{
g1->SetMinimum(20);
g1->SetMaximum(40);
}
//g1->GetYaxis()->SetLimits(-50,50);
//g1->SetMarkerStyle(21);

auto canvas1 = new TCanvas("canvas1", "Bias and sigma",1000, 500);
auto pad = new TPad("padd", "pad",0,0,1,1);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
pad->Divide(2,1,0.01,0.01);
pad->Draw();
pad->cd(1);
g->Draw("ap");
g->GetXaxis()->SetRangeUser(-0.05, 0.95);
pad->cd(2);
g1->Draw("ap");
g1->GetXaxis()->SetRangeUser(-0.05, 0.95);

if(CosmicFixed == 1){
TString nameFile1 = TString::Format("Nmix_and_Sigma_bkFixed(%d,%d,%d).pdf", static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
canvas1->SaveAs(percorso + nameFile1);
}
else{
TString nameFile1 = TString::Format("Nmix_and_Sigma(%d,%d,%d).pdf", static_cast<int>(a*100),static_cast<int>(b*100),static_cast<int>(c*100));
canvas1->SaveAs(percorso + nameFile1);
}

}
