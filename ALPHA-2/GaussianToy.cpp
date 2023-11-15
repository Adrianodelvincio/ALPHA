#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
using namespace RooFit;

void GaussianToy(){



// Create observables
RooRealVar x("x", "x", -5, 15);
double media0 = 1, media1 = 1+1.41/31.6 , varianza0 = 1, varianza1 = 1 ;
double peso0 = 0.5, peso1 = 0.5;  int Nsample = 500, Nloop = 1000;

TH1* h0 = new TH1I("h0", "N0 - N0true",50,-200,200);
TH1* h1 = new TH1I("h1", "N1 - N1true",50,-200,200); 

// Create Gaussian
RooRealVar sigma0("sigma", "sigma", varianza0,varianza0,varianza0);
RooRealVar mean0("mean", "mean", media0,media0,media0);
RooRealVar sigma1("sigma1","sigma1",varianza1,varianza1,varianza1);
RooRealVar mean1("mean1", "mean1",media1,media1,media1);
RooGaussian gaussiana("gauss0", "g0",x, mean0, sigma0);
RooGaussian gaussiana1("gauss1","g1",x, mean1, sigma1);


RooRealVar N0("N0", "N0",peso0*Nsample,0, Nsample);
RooRealVar N1("N1", "N1",peso1*Nsample,0, Nsample);
RooRealVar N00("N00", "N00",0,Nsample);
RooRealVar N11("N11", "N11",0,Nsample);

//Creo distribuzione totale, somma di due gaussiane
RooAddPdf modello("model", "model", RooArgList{gaussiana,gaussiana1}, RooArgList{N0,N1});
RooPlot *frame1 = x.frame(Title("Sum Of Pdfs"));
RooAddPdf modellofit("modelfit", "modelfit", RooArgList{gaussiana,gaussiana1}, RooArgList{N00,N11});
modello.plotOn(frame1, LineColor(28));
modello.paramOn(frame1);



double fit0 = 0,fit1 = 0,fit5 = 0;

//RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
//LOOP ON THE FIT NOW
for(int i = 0; i < Nloop; i++){
	
	// Generate a sample of Nsample events with sigma=3
	std::unique_ptr<RooDataSet> dataLoop{modello.generate(x, Nsample)};
	modellofit.fitTo(*dataLoop, Save(), PrintLevel(-1));
	//save and store the fit values
	fit0 = (N00.getVal() - Nsample*peso0) /*/Nmix_t.errorVar()->getVal()*/;
	fit1 = (N11.getVal() - Nsample*peso1)  /*/Nuw_a.errorVar()->getVal() */;
	//fit3 = N1.getVal();
	//fit4 = N1.getVal();
	//Get the ChiSquare
	RooPlot *frameLoop = x.frame(Title("data"));
	dataLoop->plotOn(frameLoop);
	modello.plotOn(frameLoop);
	//fill histogram
	h0->Fill(fit0);
	h1->Fill(fit1);
	//hchisq->Fill((frameLoop->chiSquare())*(30));
	if (i%10 == 0){
	std::cout << "Event LOOP N°" << i << std::endl;
	std::cout << "fit:      " << "N0: " << N00.getVal() << " N1: " << N11.getVal() << " Nbk " << std::endl;
	std::cout << "Expected: " << "N0: " << Nsample*peso0 << " N1: " << Nsample*peso1 << std::endl;
	std::cout << "fit - True: " << "N0: " << fit0 << " N1: " << fit1 << std::endl;}
}

//Generate a dataset of N events in x from model_analytic pdf
std::unique_ptr<RooDataSet> data{modello.generate(x,Nsample)}; //generate the data
RooDataHist *histXgen = data->binnedClone(); //create a binned histogram
modellofit.fitTo(*data, Save());
RooPlot *frame2 = x.frame(Title("Toy Model"));
histXgen->plotOn(frame2);
modellofit.plotOn(frame2);
modellofit.paramOn(frame2);

auto canvas1 = new TCanvas("d1","d1",800,800);
frame1->Draw();

auto canvas2 = new TCanvas("d2", "Toy Model data",800,800);
frame2->Draw();

auto canvas4 = new TCanvas("d4", "Toy Result", 800,800);
auto pad = new TPad("pad", "pad",0,0,1,1);
pad->Divide(2,1,0.02,0.02);
pad->Draw();

pad->cd(1);
h0->SetLineColor(1);
h0->SetLineWidth(2);
h0->Draw();
pad->cd(2);
h1->SetLineColor(1);
h1->SetLineWidth(2);
h1->Draw();
}
