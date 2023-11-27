#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
using namespace RooFit;

double algorithm(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);
double doubleThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);

void AnalysisLineShape(int start = 0, int stop = 999){

	gInterpreter->GenerateDictionary("ReadFiles", "../Headers/AnalysisLineShape.h");
	std::vector<std::string> FileList;
	TString directory = "Quadratic/";
	FileList = getFiles(start,stop, directory);
	//std::cout << FileList[0] << std::endl;
	ROOT::RDataFrame rdf("myTree", FileList);
	auto histF = rdf.Filter("Type == 0").Filter("frequence <= 80").Histo1D({"Counts","Frequence",60u,-40,60}, "frequence");
	auto histF2 = rdf.Filter("Type == 0").Filter("frequence >= 80").Histo1D({"Counts","Frequence",60u,80,200}, "frequence");
	auto histF3 = rdf.Filter("frequence <= 80").Histo1D({"Counts","Pdf1",60u,-40,60}, "frequence");
	auto histF4 = rdf.Filter("frequence >= 80").Histo1D({"Counts","Pdf2",60u,80,200}, "frequence");

	double CosmicMean = 2.0;
	double threshold = 5*CosmicMean;	// threshold considering the cosmic background
	vector<double> onset1v,onset2v;
	// IMPLEMENTING THE TOY FOR THE ALGORITHM
	for(int i = 0; i < FileList.size(); i += 2){
		std::cout << "Analizzo DataFrame " << i << std::endl;
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});		// Load i-th dataset
		auto Spectra1 = frame.Filter("frequence <= 80").Histo1D({"Counts","Frequence", 60u,-40, 60}, "frequence");
		auto Spectra2 = frame.Filter("frequence >= 80").Histo1D({"Counts","Frequence", 60u, 80, 160}, "frequence");
		auto d = frame.Display({"frequence","random", "Type", "radius"},1); d->Print();
				
		double onset1;		// Reconstructed onset and bin of the onset
		double onset2;
		
		onset1 = algorithm(Spectra1,threshold);
		onset2 = algorithm(Spectra2,threshold);
		onset1v.push_back(onset1);
		onset2v.push_back(onset2);
	}
	//histF3->Scale(1./histF3->Integral(), "width");
	//histF4->Scale(1./histF4->Integral(), "width");
	
	histF3->Scale(1./(FileList.size()/2), "width");
	histF4->Scale(1./(FileList.size()/2), "width");
	auto b = new TCanvas("b1", "Counts versus Frequencies");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	histF3->DrawClone();
	pad->cd(2);
	histF4->DrawClone();
	
	
	// Create the histograms
	std::vector<double> w(onset1v.size(),1); // weights vector
	auto h1 = new TH1D("h1","onset_{algorithm} - onset_{true}",30,-15,15);
   	h1->FillN(onset1v.size(),onset1v.data(), w.data());
   	auto h2 = new TH1D("h2","onset_{algorithm} - onset_{true}",30,110,140);
   	h2->FillN(onset2v.size(),onset2v.data(), w.data());
   	
   	auto canvas = new TCanvas("d", "Toy Result", 1000,550);
	auto pad2 = new TPad("pad2", "pad2",0,0,1,1);
	gStyle->SetOptStat(1);
	//gStyle->SetOptFit(1);
	pad2->Divide(2,1,0.005,0.005);
	pad2->Draw();
	pad2->cd(1);
	h1->SetLineColor(1);
	h1->SetLineWidth(2);
	h1->Draw();
	//legend->Draw();
	pad2->cd(2);
	h2->SetLineColor(1);
	h2->SetLineWidth(2);
	h2->Draw();
	//legend->Draw();
	canvas->SaveAs("OnsetResults.pdf");


}


double algorithm(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		if(histpdf->GetBinContent(i) >= threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
	}
	std::cout << "Threshold: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double doubleThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) >= threshold){
			if(histpdf->GetBinContent(i+1) >= threshold){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	std::cout << "		i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}
