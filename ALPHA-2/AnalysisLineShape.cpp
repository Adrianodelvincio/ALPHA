#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
using namespace RooFit;


void AnalysisLineShape(double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5){
	
	double pMix_c = Mix_c;
	double pMix_d = Mix_d;
	double c = C;	
	TString filename = TString::Format("LineShape/ToyShape1_%d_%d_c%d.root",
	static_cast<int>(pMix_c*100),
	static_cast<int>(pMix_d*100),
	static_cast<int>(c*100));
	TString filename2 = TString::Format("LineShape/ToyShape2_%d_%d_c%d.root",
	static_cast<int>(pMix_c*100),
	static_cast<int>(pMix_d*100),
	static_cast<int>(c*100));
	
	if(gSystem->AccessPathName(filename)){
		TString execute = TString::Format(".x toyLineShape.cpp(%f,%f,%f)",
		static_cast<double>(pMix_c),
		static_cast<double>(pMix_d),
		static_cast<double>(c));
		gROOT->ProcessLine(execute);
	}
	
	std::string file1(filename.Data());
	std::string file2(filename2.Data());
	ROOT::RDataFrame rdf("myTree", {file1,file2});
	auto histoF = rdf.Filter("Type == 0").Filter("frequence <= 3.6").Histo1D({"Counts","Frequence",90u,-1.2,3.6}, "frequence");
	auto histoF2 = rdf.Filter("Type == 0").Filter("frequence >= 3.6").Histo1D({"Counts","Frequence",90u,1419.2,1424.0}, "frequence");
	
	auto d1 = rdf.Display({"id", "frequence", "Type", "radius"}, 10); d1->Print();
	
	auto b = new TCanvas("b1", "Counts versus Frequencies");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	histoF->DrawClone();
	pad->cd(2);
	histoF2->DrawClone();
}
