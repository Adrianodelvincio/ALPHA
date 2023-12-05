#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"

using namespace RooFit;

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);
double firstWithVeto(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);
double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf);
double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf);

void LineShapeAnalysis(TString directory = "linear/",
					int start = 0,
					int stop = 999,
					TString ConfFile = "ToyConfiguration.txt",
					double mu = 3){
	
	gInterpreter->GenerateDictionary("ToyParser","../Headers/ConfigurationParser.h");
	gInterpreter->GenerateDictionary("ReadFiles", "../Headers/AnalysisLineShape.h");
	
	ReadConfFile Params(ConfFile);
	Params.Print();
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
	double FrequencyStep = Params.FrequencyStep;
	double startPdf1 = Params.x_cb_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.x_da_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep d-a
	double SweepStep = Params.SweepStep;
	
	std::vector<std::string> FileList;
	FileList = getFiles(start,stop, directory);
	//std::cout << FileList[0] << std::endl;
	ROOT::RDataFrame rdf("myTree", FileList);
	auto histF3 = rdf.Filter("runNumber == 1").Filter("frequence <= 1000").Histo1D({"Counts"," c to b",static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
	auto histF4 = rdf.Filter("runNumber == 1").Filter("frequence >= 1000").Histo1D({"Counts"," d to a",static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");


	double threshold = mu*CosmicBackground;	// threshold considering the cosmic background
	vector<double> onset1v,onset2v;
	vector<double> deltaOnset;
	
	// IMPLEMENTING THE TOY FOR THE ALGORITHM
	for(int i = 0; i < FileList.size(); i += 2){
		std::cout << "Analizzo DataFrame " << i << std::endl;
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});		// Load i-th dataset
		auto Spectra1 = frame.Filter("runNumber == 1")
							 .Filter("frequence <= 1000")
							 .Filter("type != 2")
							 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
		
		auto Spectra2 = frame.Filter("runNumber == 1")
							 .Filter("frequence >= 1000")
							 .Filter("type != 2")
							 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");
		
		//auto d = frame.Display({"frequence","random", "type", "radius", "delay"},1); d->Print();
		auto frame1 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Mean<double>("delay");
		auto frame2 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Mean<double>("delay");
		double delay1 = frame1.GetValue();
		double delay2 = frame2.GetValue();
		
		std::cout << "delay1: " << delay1 << " delay2: " << delay2 << std::endl;
		double onset1;			// Reconstructed onset and bin of the onset
		double onset2;
		onset1 = algorithm_2017(Spectra1);
		onset2 = algorithm_2017(Spectra2);
		//onset1 = reverse_2017(Spectra1);
		//onset2 = reverse_2017(Spectra2);
		onset1v.push_back(onset1 - (Params.x_cb_start));
		onset2v.push_back(onset2 - (Params.x_da_start));
		deltaOnset.push_back(onset2 - onset1 - (Params.x_da_start - Params.x_cb_start));
	}

	//histF3->Scale(1./histF3->Integral(), "width");
	//histF4->Scale(1./histF4->Integral(), "width");
	histF3->Scale(1./(stop));
	histF4->Scale(1./(stop));

	auto b = new TCanvas("b1", "Spectral lines");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	histF3->SetMarkerStyle(21);
	histF3->SetMarkerColor(2);
	histF3->SetLineColor(4);
	histF3->DrawClone();
	pad->cd(2);
	histF4->SetMarkerStyle(21);
	histF4->SetMarkerColor(2);
	histF4->SetLineColor(4);
	histF4->DrawClone();
	
	
	// Create the histograms
	std::vector<double> w(onset1v.size(),1); // weights vector
	auto h1 = new TH1D("h1","onset_{algorithm} - onset_{true}", 101, -50.5, 50.5);
   	h1->FillN(onset1v.size(),onset1v.data(), w.data());
   	auto h2 = new TH1D("h2","onset_{algorithm} - onset_{true}", 101, -50.5 , 50.5);
   	h2->FillN(onset2v.size(),onset2v.data(), w.data());
   	auto h3 = new TH1D("h2","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",101, -50.5 , 50.5);
   	h3->FillN(deltaOnset.size(), deltaOnset.data(), w.data());
   	
   	auto canvas = new TCanvas("d", "Toy Result", 1000,550);
	auto pad2 = new TPad("pad2", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	//gStyle->SetOptFit(1);
	pad2->Divide(2,2,0.000001,0.000001);
	pad2->Draw();
	pad2->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
	//legend->Draw();
	pad2->cd(2);
	h2->GetXaxis()->SetTitle("frequency [kHz]");
	h2->SetLineColor(2);
	h2->SetLineWidth(2);
	h2->Draw();
	pad2->cd(3);
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	//legend->Draw();
	TString name = "onsetResult"; TString endname = ".pdf"; 
	int numero = 0;
	TString folder = "Plot/";
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		name = TString::Format("onsetResult%d", numero);
	}
	canvas->SaveAs(folder + name + endname);
	
	
	/*
	auto canvas1 = new TCanvas("d1", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	*/
}

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
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

double firstWithVeto(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
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

double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) > 0){
			if(histpdf->GetBinContent(i+1) > 1){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	std::cout << "		i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){
	// bin 0 is underflow
	double onset = 0;			// onset value 
	double bin = histpdf->GetNbinsX();	// bin onset
	for(int i = histpdf->GetNbinsX(); i > 1; --i){
		if(histpdf->GetBinContent(i) < 2){
			if(histpdf->GetBinContent(i-1) < 1){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	std::cout << "		i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

