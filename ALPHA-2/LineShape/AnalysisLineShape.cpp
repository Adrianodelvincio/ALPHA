#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
using namespace RooFit;


void AnalysisLineShape(double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5, bool Execute = false){
	double pMix_c = Mix_c;
	double pMix_d = Mix_d;
	double c = C;
	if(Execute){
		TString filename = TString::Format("ToyShape1_%d_%d_c%d.root",
		static_cast<int>(pMix_c*100),
		static_cast<int>(pMix_d*100),
		static_cast<int>(c*100));
		TString filename2 = TString::Format("ToyShape2_%d_%d_c%d.root",
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
	}
	TString filename = TString::Format("LoopData_1.root");
	TString filename2 = TString::Format("LoopData_2.root");
	std::string file1(filename.Data());
	std::string file2(filename2.Data());
	ROOT::RDataFrame rdf("myTree", {file1,file2});
	auto histF = rdf.Filter("Type == 0").Filter("frequence <= 3.6").Histo1D({"Counts","Frequence",120u,-1.2,3.6}, "frequence");
	auto histF2 = rdf.Filter("Type == 0").Filter("frequence >= 3.6").Histo1D({"Counts","Frequence",120u,1419.2,1424.0}, "frequence");
	auto histF3 = rdf.Filter("frequence <= 3.6").Histo1D({"Counts","Pdf1",60u,-1.2,2}, "frequence");
	auto histF4 = rdf.Filter("frequence >= 3.6").Histo1D({"Counts","Pdf2",60u,1419.2,1424.0}, "frequence");
	auto d1 = rdf.Display({"id", "frequence", "Type", "radius"}, 10); d1->Print();

	//histF3->Scale(1./histF3->Integral(), "width");
	//histF4->Scale(1./histF4->Integral(), "width");
	auto b = new TCanvas("b1", "Counts versus Frequencies");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	histF3->DrawClone();
	pad->cd(2);
	histF4->DrawClone();
}
