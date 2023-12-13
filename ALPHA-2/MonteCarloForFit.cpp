#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
#include "TGraphErrors.h"
#include "Headers/MontecarloForFit.h"

using namespace RooFit;

void MonteCarloForFit( TString ConfFile = "configToyModel.txt"){

ReadConfFile conf(ConfFile);	// Rad the configuartion file of the program
conf.Print();					// Print the Parameters of the Simulation
gRandom->SetSeed(6); 			// Define a seed

//CREATE AN ISTOGRAM TO SAVE THE QUANTITIES OF THE MONTECARLO
TH1* hmix = new TH1I("h1", "Annihilation on trap walls : (N_{fit} - N_{gen})/#sigma",30, -6,6);
TH1* huw =  new TH1I("h2", "Residual Gas : (N_{fit} - N_{gen})/#sigma",30,-6,6);
TH1* hbk =  new TH1I("h3", "Cosmic Background : (N_{fit} - N_{gen})/#sigma",30,-6,6);
TH1* hchisq = new TH1I("hchisq", "chiquadro", 30, 0, 50);

RooRealVar x("x","r [cm]",0.,4.);
x.setBins(30);

// ANNIHILATION ON WALL PDF
RooRealVar mu("mu", "mu", conf.mu);
RooRealVar sigMix("sigMix", "sigMix", conf.sigWall);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

//PDF UWlosses
RooRealVar sigRay("sigRay", "sigma", conf.sigRay);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

//PDF Cosmic
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));

RooRealVar Nmix_t("Nmix","Nmix", 	conf.N * conf.a, -3000, +3000);
RooRealVar Nuw_a ("Ngas", "Ngas",	conf.N * conf.b, -3000, +3000);
RooRealVar Nbk_a ("Ncosmic", "Ncosmic", conf.Ncosmic, -3000, +3000);
RooPlot *frame1 = x.frame(Title("Sum Of Three Pdfs"));
//Model to generate the data
RooAddPdf model_gen("model", "model", RooArgList{gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});
model_gen.plotOn(frame1, LineColor(28));
model_gen.paramOn(frame1);

//Generate a dataset of N events in x from model_analytic pdf
std::unique_ptr<RooDataSet> data{model_gen.generate(x,Extended())}; //generate the data
RooDataHist *histXgen = data->binnedClone(); //create a binned histogram

// FIT the toy model
RooRealVar Nmix_f("Nfit_{mix}","Nmix",	conf.a*conf.N,-3000,+3000);
RooRealVar Nuw_f("Nfit_{gas}","Ngas",	conf.b*conf.N,-3000,+3000);
RooRealVar Nbk_f("Nfit_{cosmic}", "Nbk",conf.Ncosmic, -3000,+3000);

RooAddPdf model_forfit("model2", "model", RooArgList{gauss_Mix,Rayleigh,linearFit},RooArgList{Nmix_f,Nuw_f,Nbk_f});
double fit0, fit1, fit2;
double fit3 = 0,fit4 = 0,fit5 = 0;
//Disable some unuseful print
RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

//LOOP ON THE FIT NOW
for(int i = 0; i < conf.Nloop; i++){
	//generate the data
	std::unique_ptr<RooDataSet> dataLoop{model_gen.generate(x,Extended())};
	model_forfit.fitTo(*dataLoop, Save(), PrintLevel(-1));
	//save and store the fit values
	fit0 = (Nmix_f.getVal() - conf.N*conf.a)/Nmix_f.errorVar()->getVal();
	fit1 = (Nuw_f.getVal()  - conf.N*conf.b)/Nuw_f.errorVar()->getVal();
	fit2 = (Nbk_f.getVal()  - conf.Ncosmic )/Nbk_f.errorVar()->getVal();
	//Get the ChiSquare
	RooPlot *frameLoop = x.frame(Title("data"));
	dataLoop->plotOn(frameLoop);
	model_forfit.plotOn(frameLoop);
	//fill histogram
	hmix->Fill(fit0);
	huw->Fill(fit1);
	hbk->Fill(fit2);
	hchisq->Fill((frameLoop->chiSquare())*(30));
	if(i%10 == 0){
		std::cout << "Event LOOP N°" << i << " Gen Event: " << dataLoop->numEntries() <<std::endl;
		std::cout << "fit:      " << "Nmix: " << Nmix_f.getVal() << " Nuw: " << Nuw_f.getVal() << " Nbk " << Nbk_f.getVal() << std::endl;
		std::cout << "Expected: " << "Nmix: " << conf.N*conf.a << " Nuw: " << conf.N*conf.b << " Nbk: " << conf.Ncosmic << std::endl;
	}
}

//Fit the Histograms of N - Ntrue
hmix->Fit("gaus");
huw->Fit("gaus");
hbk->Fit("gaus");

auto legend = new TLegend(0.1,0.7,0.48,0.9);
legend->SetHeader("PDF Composition","C"); // option "C" allows to center the header
TString coeffMix = TString::Format("events on walls: %d", static_cast<int>(conf.a*conf.N));
TString coeffUw = TString::Format("events res. gas: = %d", static_cast<int>(conf.b*conf.N));
TString coeffbk = TString::Format("cosmics = %.1f", conf.Ncosmic);
legend->AddEntry(hmix,coeffMix);
legend->AddEntry(hmix,coeffUw);
legend->AddEntry(hmix,coeffbk);
legend->SetTextSize(.026);

// LAYOUT OF THE PLOTS
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
gStyle->SetFitFormat("2.3g");
// hmix layout
hmix->SetTitleSize(0.1,"t");
hmix->SetLineColor(1);
hmix->SetLineWidth(2);
hmix->SetAxisRange(0.,hmix->GetMaximum() + 40,"Y");
auto function0 = hmix->GetFunction("gaus");
function0->SetLineStyle(9);
// huw layout
huw->SetTitleSize(0.1,"t");
huw->SetLineColor(1);
huw->SetLineWidth(2);
huw->SetAxisRange(0., huw->GetMaximum() + 40,"Y");
auto function1 = huw->GetFunction("gaus");
function1->SetLineStyle(9);
// hbk layout
hbk->SetTitleSize(0.95,"t");
hbk->SetLineColor(1);
hbk->SetLineWidth(2);
hbk->SetAxisRange(0., hbk->GetMaximum() + 40,"Y");
auto function2 = hbk->GetFunction("gaus");
function2->SetLineStyle(9);

// Plot the Result
auto canvas4 = new TCanvas("d4", "Toy Result", 1500,800);
auto pad = new TPad("pad", "pad",0,0,1,1);

pad->Divide(3,1,0.01,0.01);
pad->Draw();
pad->cd(1);
hmix->Draw();
legend->Draw();
pad->cd(2);
huw->Draw();
pad->cd(3);
hbk->SetTitleSize(0.1,"t");
hbk->Draw();

TString percorso = TString::Format("PlotMLEfit/N%d/", conf.N);
TString nameFile0 = TString::Format("Reconstructed_(%d,%d,%.1f).pdf",static_cast<int>(conf.a*conf.N),static_cast<int>(conf.b*conf.N),conf.Ncosmic);
canvas4->SaveAs(percorso + nameFile0);

auto canvas5 = new TCanvas("d5","Chi square distribution",800,800);
hchisq->SetLineWidth(2);
hchisq->Draw();


auto canvas6 = new TCanvas("d6", "Toy Result", 1000,700);
auto pad1 = new TPad("pad", "pad",0,0,1,1);
gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
pad1->Divide(2,1,0.01,0.01);
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

TString nameFile1 = TString::Format("Reconstruced_(%d,%d).pdf",static_cast<int>(conf.a*conf.N),static_cast<int>(conf.b*conf.N));
canvas6->SaveAs(percorso + nameFile1);

}


void SecondMontecarlo(){
double a = 0.9;  // weight mix
double b = 0. ;  // weight residual gas
double c = 0.10; // weight cosmic 
int CosmicFixed = 0;
int N = 1000;    // Number of event to be generated
int Npoint = 30; // Number of step in weight scan
int Nnested = 100;// number of trials of the fit

//Generate dictionary
gInterpreter->GenerateDictionary("MontecarloForFit", "Headers/MontecarloForFit.h");

// Fix the cosmic Background
double Nfix = 10.2;
RooRealVar x("x","r [cm]",0.,4.); x.setBins(30); // Set the variable to store the data

RooRealVar mu("mu", "mu", 2.38,2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85,0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

RooRealVar sigRay("sigRay", "#sigma_{rayleigh}", 1.722,1.722); // fixing sigma from the fit value
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x)); //PDF Cosmic

//Model to generate the data
RooRealVar Nmix_t("Nmix",    "Nmix", 	1, -(2*N), +(2*N));
RooRealVar Nuw_a ("Ngas",    "Ngas", 	1, -(2*N), +(2*N));
RooRealVar Nbk_a ("Ncosmic", "Ncosmic", 1, -(2*N), +(2*N));
RooAddPdf model_gen("model", "model", RooArgList{ gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix_t,Nuw_a,Nbk_a});

//Model for the Fit
RooRealVar Nmix_f("Nmix fit",    "Nmix",    1,-(2*N), +(2*N));
RooRealVar Nuw_f ("Ngas fit",    "Ngas",    1,-(2*N), +(2*N));
RooRealVar Nbk_f ("Ncosmic fit", "Ncosmic", 1,-(2*N), +(2*N));
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
	double wmix = 0 + static_cast<double>(i)*(1 - c)/(Npoint-1); //Compute the weight wmix
	weight.push_back(wmix); // store the weight
	ChangeWeight(&Nmix_t, &Nuw_a, &Nbk_a, wmix, c, N); // Update the expected events
	for(int j = 0; j < Nnested; j++){
		std::unique_ptr<RooDataSet> dataLoop{model_gen.generate(x,Extended())}; // generate data
		model_forfit.fitTo(*dataLoop, Save(), PrintLevel(-1)); 			// fit data
		SetVectors(j,mix_l,uw_l,bk_l,Nmix_f,Nuw_f,Nbk_f,Nmix_t,Nuw_a,Nbk_a); 	// store parameters
		mixErrors_l.push_back(Nmix_f.errorVar()->getVal()); 			//set error of the parameters
	if(j%10 == 0 ){std::cout << "Event gen " << dataLoop->numEntries() << std::endl;}
	}
	//save and store the fit values
	mix.push_back(media(mix_l));
	uw.push_back(media(uw_l));
	bk.push_back(media(bk_l));
	mixErrors.push_back(media(mixErrors_l));
	Errors.push_back(media(mixErrors_l)/sqrt(Nnested));
	mix_l.clear(); uw_l.clear(); bk_l.clear(); mixErrors_l.clear(); // Clear the vectors
	if(i%1 == 0){PrintInfo(Nmix_f,Nuw_f,Nbk_f,Nmix_t, Nuw_a, Nbk_a, wmix, i);} // Print some info
}

//Visualize the Montecarlo
vector<double> errorx(30,0); // define 0 error on x axis
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
