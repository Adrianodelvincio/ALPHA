#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include <TMath.h>
using namespace RooFit;

void AnalysisMLE(){

// DATA TO FIT
TString cartella = TString::Format("Spectroscopy/Dataset/");
TString filename = TString::Format("r68465_uw_exp_freq4.vertex.csv");

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
mix2_rdf.Snapshot("myTree", "Spectroscopy/RootCut1Data/MixCut1.root", {"X","Y","Radius"});

RooRealVar x("x", "r [cm]", 0, 4);
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
uw2_rdf.Snapshot("myTree", "Spectroscopy/RootCut1Data/UWCut1.root", {"X","Y","Radius"});

RooDataHist Uw_h("dh1", "dh1", x, Import(*histUw));
RooPlot *frame2 = x.frame(Title("UW losses PDF"));
frame2->GetYaxis()->SetTitle("Counts");
// EXTRACT THE PDF FROM THE HISTOGRAM
RooHistPdf PdfUwlosses("uwlossespdf", "uwlossespdf", x, Uw_h, 0);
Uw_h.plotOn(frame2);
PdfUwlosses.plotOn(frame2);

// FIT THE DATA WITH A RAYLEIGH

RooRealVar sigRay("sigRay", "sigma", 1.722);
RooRealVar sigRay2("sigRay", "sigma", 1.722,0,100);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay)))", RooArgSet(x,sigRay));
//RUN THE FIRST TIME FOR THE PARAMETERS OF THE RAYLEIGH

RooGenericPdf Rayleigh2("line", "linear model", " TMath::Abs(x/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay)))", RooArgSet(x,sigRay2));
Rayleigh2.fitTo(Uw_h); // FIT
Rayleigh2.plotOn(frame2, LineColor(kRed));
Double_t chi2ray = frame2->chiSquare();
std::cout << "chi square rayleigh: " << frame2->chiSquare() << std::endl; 
TString chis2ray = TString::Format("Chisquare = %.1f ndof %d", chi2ray*27,27);
Rayleigh2.paramOn(frame2,Label(chis2ray));

auto canvas1 = new TCanvas("d1", "d1",800,800);
frame2->Draw();

//COSMIC

ROOT::RDataFrame cosmic_rdf("myTree",{"DataSetROOT/r68949_cosmics.vertex.root",
	"DataSetROOT/r69177_cosmics.vertex.root",
	"DataSetROOT/r69207_cosmics.vertex.root",
	"DataSetROOT/r69219_cosmics.vertexroot"});

auto cosmic2_rdf = cosmic_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto histCosmic = cosmic2_rdf.Histo1D({"Cosmic","Counts",30u,0.,4.}, "Radius");
cosmic2_rdf.Snapshot("myTree", "Spectroscopy/RootCut1Data/CosmicCut1.root", {"X","Y","Radius"});

RooDataHist cosm_h("dh2", "dh2", x, Import(*histCosmic));
RooPlot *frame3 = x.frame(Title("Cosmic PDF"));

//EXTRACT THE PDF FROM THE PLOT
RooHistPdf PdfBk("Cosmic", "cosmic", x, cosm_h, 0);
cosm_h.plotOn(frame3);
PdfBk.plotOn(frame3);

// FIT THE PDF WITH A LINE

RooRealVar m("m", "m", 0,1000000);
//RooRealVar q("q", "q",0.25, 0.25 - 200);
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));
RooAddPdf line_pdf("modell","modell", RooArgList{linearFit}, RooArgList{m} );
line_pdf.fitTo(cosm_h, PrintLevel(-1));
line_pdf.paramOn(frame3);
line_pdf.plotOn(frame3, LineColor(kRed));

auto canvas3 = new TCanvas("d2","d2",800,800);
frame3->Draw();


// PRINT ALL THE PDF TOGHETHER IN THE SAME CANVAS

auto canvasPdf = new TCanvas("p0", "PDF normalized", 800, 800);
RooPlot *allpdf = x.frame(Title("All Pdf"));

PdfBk.plotOn(allpdf, LineColor(kBlue), Name("Cosmic"));
PdfMixing.plotOn(allpdf, LineColor(kRed), Name("Mixing"));
PdfUwlosses.plotOn(allpdf, LineColor(kGreen), Name("Uw"));

// Aggiungo una legenda al plot
TLegend *leg1 = new TLegend(0.65,0.73,0.86,0.87);
leg1->SetFillColor(kWhite);
leg1->SetLineColor(kWhite);
leg1->AddEntry(allpdf->findObject("Cosmic"), "Cosmic Cemplate", "LP");
leg1->AddEntry(allpdf->findObject("Mixing"), "Mixing Template", "LP");
leg1->AddEntry(allpdf->findObject("Uw"), "Uw Template","LP");
allpdf->Draw();
leg1->Draw();

auto f4_rdf = ROOT::RDF::MakeCsvDataFrame(cartella + filename);
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

// Fit with the sum of all distributions
RooRealVar Nmix("Nmix", "Nmix", 0., 200);
RooRealVar Nuw("Nuw", "Nuw", -100, 200);
RooRealVar Nbk("Nbk", "Nbk", -100, 200);
RooPlot *frame4 = x.frame(Title("Fit to f4 Dataset"));
RooAddPdf model("model","model", RooArgList{PdfMixing,PdfUwlosses,PdfBk}, RooArgList{Nmix, Nuw, Nbk} );
model.fitTo(f4_h, Extended(),PrintLevel(-1));
f4_h.plotOn(frame4);
model.plotOn(frame4, LineColor(kRed));

// Print the chisquare on the plot
Double_t chi2 = frame4->chiSquare();
TString chis2Line = TString::Format("Chisquare = %f ", chi2);
model.paramOn(frame4, Label(chis2Line));

auto canvas4 = new TCanvas("d4","d4", 800, 800);
// Calculate the Chisquare
std::cout << "chi^2 = " << frame4->chiSquare() << std::endl;
frame4->Draw();

// FIT WITH COSMIC FIXED
RooRealVar Nmix_f("Nmix", "Nmix", 0., 200);
RooRealVar Nuw_f("Nuw", "Nuw", -100, 100);
RooRealVar Nbk_f("Nbk", "Nbk", 10,10);
RooPlot *frame5 = x.frame(Title("Fit, Cosmic Fixed"));
RooAddPdf model_bkfixed("model","model", RooArgList{PdfMixing,PdfUwlosses,PdfBk}, RooArgList{Nmix_f, Nuw_f, Nbk_f} );
model_bkfixed.fitTo(f4_h,PrintLevel(-1)); // FIT
f4_h.plotOn(frame5);
model_bkfixed.plotOn(frame5, LineColor(kBlue));

// Print the chisquare on the plot
Double_t chi2_bkfixed = frame5->chiSquare();
TString chis2Line_fixed = TString::Format("Chisquare = %f ", chi2_bkfixed);
model_bkfixed.paramOn(frame5, Label(chis2Line_fixed));

auto canvas6 = new TCanvas("d6"," Cosmic Fixed", 800, 800);
// Calculate the Chisquare
std::cout << "chi^2 bk fixed = " << frame5->chiSquare() << std::endl;
frame5->Draw();

// FIT WITH ONLY MIXING + COSMIC(FIXED)

RooRealVar Nmix_f2("Nmix", "Nmix", 0., 200);
RooPlot *frame6 = x.frame(Title("Fit, MIxing only"));
RooAddPdf model_onlyMix("model","model", RooArgList{PdfMixing,PdfBk}, RooArgList{Nmix_f2, Nbk_f} );
model_onlyMix.fitTo(f4_h,PrintLevel(-1)); // FIT
f4_h.plotOn(frame6);
model_onlyMix.plotOn(frame6, LineColor(kGreen));

// Print the chisquare on the plot
Double_t chi2_onlyMix = frame6->chiSquare();
TString chis2Line_onlyMix = TString::Format("Chisquare = %f ", chi2_onlyMix);
model_onlyMix.paramOn(frame6, Label(chis2Line_onlyMix));

auto canvas7 = new TCanvas("d7"," Only Mix", 800, 800);
// Calculate the Chisquare
std::cout << "chi^2 only mix = " << frame6->chiSquare() << std::endl;
frame6->Draw();

// FIT WITH THE ANALYTIC FORMULA

//FOR MIXING USE A GAUSSIAN
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

RooRealVar Nmix_t("Nmix","Nmix",0.,200.);
RooRealVar Nuw_a("Nuw","Nuw", -100,100);
RooRealVar Nbk_a("Nbk", "Nbk",10.2,10.2);

RooPlot *analyticframe = x.frame(Title(filename));
analyticframe->GetYaxis()->SetTitle("Counts");
RooAddPdf model_analytic("model", "model", RooArgList{/*PdfMixing*/ gauss_Mix,Rayleigh,linearFit}, RooArgList{Nmix,Nuw_a,Nbk_a});
model_analytic.fitTo(f4_h, PrintLevel(-1));
f4_h.plotOn(analyticframe);
model_analytic.plotOn(analyticframe, LineColor(28));
TString chiquadrato = TString::Format("chisq = %.1f, ndof = %d",analyticframe->chiSquare()*(30), (30-3));
model_analytic.paramOn(analyticframe, Label(chiquadrato));

auto canvas8 = new TCanvas("d8","d8",800,800);
analyticframe->Draw();
}
