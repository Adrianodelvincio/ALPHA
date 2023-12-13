#include <iostream>
#include <TMath.h>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
using namespace RooFit;

void TemplateMLE(){
// DATA TO FIT
TString cartella = TString::Format("DataSetROOT/");
TString filename = TString::Format("r68465_f4");
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
mix2_rdf.Snapshot("myTree", "DataSetROOT/MixCut1.root", {"X","Y","Radius"});

RooRealVar x("x", "r [cm]", 0, 4);
x.setBins(30);

RooDataHist mix_h("dh0", "dh0", x, Import(*histMix));
RooPlot *frame1 = x.frame(Title("Mixing PDF"));
RooHistPdf PdfMixing("mixingpdf", "mixingpdf", x, mix_h, 0);

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
uw2_rdf.Snapshot("myTree", "DataSetROOT/UWCut1.root", {"X","Y","Radius"});

RooDataHist Uw_h("dh1", "dh1", x, Import(*histUw));
RooPlot *frame2 = x.frame(Title("GAS losses PDF"));
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
TString chis2ray = TString::Format("#chi^{2} = %.1f ndof %d", chi2ray*27,27);
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
cosmic2_rdf.Snapshot("myTree", "DataSetROOT/CosmicCut1.root", {"X","Y","Radius"});

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

// PRINT ALL THE NORMALIZED DISTRIBUTIONS TOGHETHER IN THE SAME CANVAS
RooPlot *allpdf = x.frame(Title("Normalized Distributions"));

PdfBk.plotOn(allpdf, LineColor(kBlue), Name("Cosmic"));
PdfMixing.plotOn(allpdf, LineColor(kRed), Name("Mixing"));
PdfUwlosses.plotOn(allpdf, LineColor(kGreen), Name("Gas"));

// Aggiungo una legenda al plot
TLegend *leg1 = new TLegend(0.65,0.73,0.86,0.87);
leg1->SetFillColor(kWhite);
leg1->SetLineColor(kWhite);
leg1->AddEntry(allpdf->findObject("Cosmic"), "Cosmic Distribution", "LP");
leg1->AddEntry(allpdf->findObject("Mixing"), "Annihi. on Walls Distribution", "LP");
leg1->AddEntry(allpdf->findObject("Gas"), "Residual Gas Annihi. Distribution","LP");
allpdf->SetYTitle("");

auto canvasPdf = new TCanvas("p0", "PDF normalized", 800, 800);
allpdf->Draw();
leg1->Draw();

histCosmic->Scale(1/histCosmic->Integral());
histUw->Scale(1/histUw->Integral());
histMix->Scale(1/histMix->Integral());
auto canvasPdf2 = new TCanvas("p1", "PDF radial density", 800,800);

for(int i = 1; i <= histCosmic->GetNbinsX(); i++){
	double radius = histCosmic->GetBinCenter(i);
	histCosmic->SetBinContent(i, histCosmic->GetBinContent(i)/(2*TMath::Pi()*radius));
	histUw->SetBinContent(i, histUw->GetBinContent(i)/(2*TMath::Pi()*radius));
	histMix->SetBinContent(i, histMix->GetBinContent(i)/(2*TMath::Pi()*radius));
}

gStyle->SetOptStat(0);
histUw->SetTitle("Radial Density Distributions");
histUw->SetXTitle("Radius r [cm]");
histUw->SetYTitle("Counts cm^{-1}");
histUw->SetLineColor(3);
histUw->SetLineWidth(3);
histUw->DrawClone("HIST");
histCosmic->SetLineColor(4);
histCosmic->SetLineWidth(3);
histCosmic->DrawClone("HISTsame");
histMix->SetLineColor(2);
histMix->SetLineWidth(3);
histMix->DrawClone("histsame");
//canvasPdf2->SaveAs("PlotMLEfit/RadialDensity.pdf");

// DATA FOR THE FIT f
TString formato = TString::Format(".root");
ROOT::RDataFrame  f4_rdf("myTree",{cartella + filename + formato});
auto displ = f4_rdf.Display({"CutsType1","CutsType2", "X", "Y", "Z"}, 5);
    displ->Print();

auto f4_2_rdf = f4_rdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)").Filter("CutsType1 ==  \" 1\"");
auto Histf4_h = f4_2_rdf.Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius"); 


RooDataHist f4_h("dh3", "dh3", x, Import(*Histf4_h));

// Fit with the sum of all distributions
RooRealVar Nmix("Nmix", "Nmix", 0., 200);
RooRealVar Nuw("Ngas", "Ngas", -100, 200);
RooRealVar Nbk("Ncosmic", "Ncosmic", -100, 200);
RooPlot *frame4 = x.frame(Title("Fit to f4 Dataset"));
RooAddPdf model("model","model", RooArgList{PdfMixing,PdfUwlosses,PdfBk}, RooArgList{Nmix, Nuw, Nbk} );
model.fitTo(f4_h, Extended(),PrintLevel(-1));
f4_h.plotOn(frame4);
model.plotOn(frame4, LineColor(kRed));

// Print the chisquare on the plot
Double_t chi2 = frame4->chiSquare();
TString chis2Line = TString::Format("#chi^{2} = %f ", chi2);
model.paramOn(frame4, Label(chis2Line));

auto canvas4 = new TCanvas("d4","d4", 800, 800);
// Calculate the Chisquare
std::cout << "chi^2 = " << frame4->chiSquare() << std::endl;
frame4->Draw();

}
