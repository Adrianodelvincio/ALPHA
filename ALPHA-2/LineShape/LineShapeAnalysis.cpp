#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/Algorithms.h"
using namespace RooFit;


void LineShapeAnalysis(TString directory = "linear/",
					int start = 0,
					int stop = 999,
					double mu = 3, 		// Threshold coefficient
					double fraction = 0.1, 	// constant fraction discrimination
					TString folder = "Plot/"
					){
	TString ConfFile = directory + "ToyConfiguration.txt";
	std::cout << ConfFile << std::endl;
	gInterpreter->GenerateDictionary("ToyParser","../Headers/ConfigurationParser.h");
	//gInterpreter->GenerateDictionary("ReadFiles", "../Headers/AnalysisLineShape.h");
	
	ReadConfFile Params(ConfFile);
	Params.Print();
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
	double FrequencyStep = Params.FrequencyStep;
	double startPdf1 = Params.x_cb_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.x_da_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep d-a
	double SweepStep = Params.SweepStep;
	
	std::vector<std::string> FileList;
	FileList = getFiles(start,stop, directory);
	ROOT::RDataFrame rdf("myTree", {FileList[0], FileList[1]});
	auto hist_ctob = rdf.Filter("runNumber == 0")
			.Filter("frequence <= 1000")
			.Filter("type != 2")
			.Histo1D({"Counts"," c to b",static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
	auto ctob_withCosmic = rdf.Filter("runNumber == 0")
			.Filter("frequence <= 1000")
			//.Filter("type != 2")
			.Histo1D({"Counts"," c to b",static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
	auto histF4 = 	rdf.Filter("runNumber == 0")
			.Filter("frequence >= 1000")
			.Filter("type != 2")
			.Histo1D({"Counts"," d to a",static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");

	vector<double> v1_2017, v2_2017;
	vector<double> v1_rev, v2_rev;
	vector<double> v1_thr, v2_thr;
	vector<double> v1_cfrac, v2_cfrac;
	vector<double> v1_neigh, v2_neigh;
	vector<double> diff_2017, diff_rev, diff_thr, diff_cfrac, diff_neigh;
	
	// IMPLEMENTING THE TOY FOR THE ALGORITHM
	for(int i = 0; i < FileList.size(); i += 2){
		std::cout << "Analizzo DataFrame " << i << "\n" << std::endl;
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});		// Load i-th dataset
		auto Spectra1 = frame.Filter("runNumber == 1")
							 .Filter("frequence <= 1000")
							 //.Filter("type != 2")
							 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
		
		auto Spectra2 = frame.Filter("runNumber == 1")
							 .Filter("frequence >= 1000")
							 //.Filter("type != 2")
							 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");
		
		//auto d = frame.Display({"frequence","random", "type", "radius", "lineShift"},1); d->Print();
		auto frame1 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Mean<double>("lineShift");
		auto frame2 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Mean<double>("lineShift");
		auto frame3 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Take<double>("lineShift");
		auto frame4 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Take<double>("lineShift");
		auto lineShiftda = frame3.GetValue(); std::cout << "LineShift c to b : " << lineShiftda[0] << std::endl;
		auto lineShiftcb = frame4.GetValue(); std::cout << "LineShift c to b : " << lineShiftcb[0] << std::endl;
		double delay1 = frame1.GetValue();
		double delay2 = frame2.GetValue();
		
		//std::cout << "delay1: " << delay1 << " delay2: " << delay2 << std::endl;
		double onset1;				// Reconstructed onset
		double onset2;				// Reconstructed onset
		double threshold = mu*CosmicBackground; // threshold considering the cosmic background
		
		// THRESHOLD
		onset1 = firstOverThreshold(Spectra1, mu);
		onset2 = firstOverThreshold(Spectra2, mu);
		v1_thr.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_thr.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_thr.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// 2017 FOWARD
		onset1 = algorithm_2017(Spectra1);
		onset2 = algorithm_2017(Spectra2);
		
		v1_2017.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_2017.push_back(onset2 - (Params.x_da_start + lineShiftda[0])); 
		diff_2017.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// REVERSED 2017 
		onset1 = reverse_2017(Spectra1);
		onset2 = reverse_2017(Spectra2);
		
		v1_rev.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_rev.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_rev.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// CONSTANT FRACTION
		onset1 = constFrac(Spectra1, fraction);
		onset2 = constFrac(Spectra2, fraction);
		
		v1_cfrac.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_cfrac.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_cfrac.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// SUM NEIGHBORS
		onset1 = sumNeighbors(Spectra1, CosmicBackground, mu);
		onset2 = sumNeighbors(Spectra2, CosmicBackground, mu);
		
		v1_neigh.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_neigh.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_neigh.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
	}

	//hist_ctob->Scale(1./hist_ctob->Integral(), "width");
	//ctob_withCosmic->Scale(1./ctob_withCosmic->Integral(), "width");
	//hist_ctob->Scale(1./(stop));
	//ctob_withCosmic->Scale(1./(stop));

	auto b = new TCanvas("b1", "Spectral lines");
	auto pad = new TPad("pad1", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	hist_ctob->SetTitle("c to b without cosmic events");
	hist_ctob->SetMarkerStyle(21);
	hist_ctob->SetMarkerColor(2);
	hist_ctob->SetLineColor(4);
	hist_ctob->DrawClone();
	hist_ctob->DrawClone("SAME P");
	pad->cd(2);
	ctob_withCosmic->SetTitle("c to b with cosmic events");
	ctob_withCosmic->SetMarkerStyle(21);
	ctob_withCosmic->SetMarkerColor(2);
	ctob_withCosmic->SetLineColor(4);
	ctob_withCosmic->DrawClone();
	ctob_withCosmic->DrawClone("SAME P");
	
	auto h1 = new TH1D("h1","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h2 = new TH1D("h2","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h3 = new TH1D("h2","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	
	//2017 FOWARD
	std::vector<double> w(v1_2017.size(),1); // weights vector
	
   	h1->FillN(v1_2017.size(),v1_2017.data(), w.data());
   	h2->FillN(v2_2017.size(),v2_2017.data(), w.data());
   	h3->FillN(diff_2017.size(), diff_2017.data(), w.data());
   	
   	auto canvas = new TCanvas("d", "2017 algorithm", 1000,550);
	auto pad2 = new TPad("pad2", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad2->Divide(2,2,0.001,0.001);
	pad2->Draw();
	pad2->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
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
	TString name = "2017_foward"; TString endname = ".pdf"; 
	int numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas->SaveAs(folder + name + endname);
	
	auto delta = new TCanvas("dd", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	name = "delta_2017_foward";
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta->SaveAs(folder + name + endname);
	h1->Reset("ICESM"); h2->Reset("ICESM"); h3->Reset("ICESM");
	// END 2017
	
	// 2017 REVERSED
   	h1->FillN(v1_rev.size(),v1_rev.data(), w.data());
   	h2->FillN(v2_rev.size(),v2_rev.data(), w.data());
   	h3->FillN(diff_rev.size(), diff_rev.data(), w.data());
	
	auto canvas1 = new TCanvas("d1", "2017 reversed", 1000,550);
	auto pad3 = new TPad("pad3", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad3->Divide(2,2,0.001,0.001);
	pad3->Draw();
	pad3->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
	pad3->cd(2);
	h2->GetXaxis()->SetTitle("frequency [kHz]");
	h2->SetLineColor(2);
	h2->SetLineWidth(2);
	h2->Draw();
	pad3->cd(3);
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	name = TString::Format("2017_reversed");
	
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas1->SaveAs(folder + name + endname);
	
	auto delta1 = new TCanvas("dd1", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	name = "delta_2017_reversed";
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta1->SaveAs(folder + name + endname);
	h1->Reset("ICESM"); h2->Reset("ICESM"); h3->Reset("ICESM");
	
	// THRESHOLD
	h1->FillN(v1_thr.size(),v1_thr.data(), w.data());
   	h2->FillN(v2_thr.size(),v2_thr.data(), w.data());
   	h3->FillN(diff_thr.size(), diff_thr.data(), w.data());
	
	auto canvas2 = new TCanvas("d2", "threshold", 1000,550);
	auto pad4 = new TPad("pad4", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad4->Divide(2,2,0.001,0.001);
	pad4->Draw();
	pad4->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
	pad4->cd(2);
	h2->GetXaxis()->SetTitle("frequency [kHz]");
	h2->SetLineColor(2);
	h2->SetLineWidth(2);
	h2->Draw();
	pad4->cd(3);
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	
	name = TString::Format("Threshold_%d", static_cast<int>(mu));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas2->SaveAs(folder + name + endname);
	
	auto delta2 = new TCanvas("dd2", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	name = TString::Format("delta_Threshold_%d", static_cast<int>(mu));
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta2->SaveAs(folder + name + endname);
	h1->Reset("ICESM"); h2->Reset("ICESM"); h3->Reset("ICESM");
	
	
	//COSTANT FRACTION
	h1->FillN(v1_cfrac.size(),v1_cfrac.data(), w.data());
   	h2->FillN(v2_cfrac.size(),v2_cfrac.data(), w.data());
   	h3->FillN(diff_cfrac.size(), diff_cfrac.data(), w.data());
	
	auto canvas3 = new TCanvas("d3", "constFract", 1000,550);
	auto pad5 = new TPad("pad5", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad5->Divide(2,2,0.001,0.001);
	pad5->Draw();
	pad5->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
	pad5->cd(2);
	h2->GetXaxis()->SetTitle("frequency [kHz]");
	h2->SetLineColor(2);
	h2->SetLineWidth(2);
	h2->Draw();
	pad5->cd(3);
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	
	name = TString::Format("constFract_%d", static_cast<int>(100*fraction));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas3->SaveAs(folder + name + endname);
	
	auto delta3 = new TCanvas("dd3", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	name = TString::Format("delta_constFract_%d", static_cast<int>(100*fraction));
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta3->SaveAs(folder + name + endname);
	h1->Reset("ICESM"); h2->Reset("ICESM"); h3->Reset("ICESM");
	
	// SUMNEIGHBORS	
	h1->FillN(v1_neigh.size(),v1_neigh.data(), w.data());
   	h2->FillN(v2_neigh.size(),v2_neigh.data(), w.data());
   	h3->FillN(diff_neigh.size(), diff_neigh.data(), w.data());
	
	auto canvas4 = new TCanvas("d4", "sumNeighbours", 1000,550);
	auto pad6 = new TPad("pad6", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad6->Divide(2,2,0.001,0.001);
	pad6->Draw();
	pad6->cd(1);
	h1->GetXaxis()->SetTitle("frequency [kHz]");
	h1->SetLineColor(2);
	h1->SetLineWidth(2);
	h1->Draw();
	pad6->cd(2);
	h2->GetXaxis()->SetTitle("frequency [kHz]");
	h2->SetLineColor(2);
	h2->SetLineWidth(2);
	h2->Draw();
	pad6->cd(3);
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	
	name = TString::Format("sumNeighbors_sigma%d", static_cast<int>(mu));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas4->SaveAs(folder + name + endname);
	
	auto delta4 = new TCanvas("dd4", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h3->GetXaxis()->SetTitle("frequency [kHz]");
	h3->SetLineWidth(2);
	h3->SetLineColor(38);
	h3->Draw();
	name = TString::Format("delta_sumNeighbors_sigma%d", static_cast<int>(mu));
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta4->SaveAs(folder + name + endname);
	h1->Reset("ICESM"); h2->Reset("ICESM"); h3->Reset("ICESM");
	//TString name = TString::Format("sumNeighbours_sigma%d(cosmic=0)", static_cast<int>(mu)); TString endname = ".pdf"; 

}
