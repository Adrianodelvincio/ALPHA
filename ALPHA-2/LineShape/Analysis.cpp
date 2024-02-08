#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/Algorithms.h"
using namespace RooFit;


void  Analysis(	TString directory,
		int start,              // analyze trials starting from start
		int stop,		// analyze trials until stop
		double Nfilter,		// Filter parameter
		double fraction,	// Parameter of Constant Fraction
		double Nsigma,		// Paratemer for t-student test
		double Nthr,		// Parameter of Threshold algorithm
		double thr1,            // first threshold for forward
		double thr2,            // second threshold for forward
		TString folder = "Plot/"
		){

	TString ConfFile = directory + "ToyConfiguration.txt";
	std::cout << ConfFile << std::endl;
	
	// Read configuration values
	ReadConfFile Params(ConfFile); // Read the values from the configuration file 
	// Print configuration values
	Params.Print();

	std::vector<std::string> FileList;
	FileList = getFiles(start,stop, directory); // file list to be analyzed
	
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Baseline
	double FrequencyStep = Params.FrequencyStep;

	// JUST FOR VISUALIZATION, uncomment to watch the data
	// Show the first two files
	
	double startPdf1 = Params.cb_start - (Params.FrequencyStep)*(6 + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.da_start - (Params.FrequencyStep)*(6 + 0.5);	// Start of frequency sweep d-a
	
	ROOT::RDataFrame rdf("myTree", {FileList[0], FileList[1]});
	auto hist_ctob = rdf.Filter("repetition == 0")
			.Filter("mwfrequence <= 1000")
			//.Filter("type != 2")
			.Histo1D({"Counts"," c to b",static_cast<int>(50),startPdf1, startPdf1 + 50*Params.FrequencyStep }, "mwfrequence");
	auto ctob_withCosmic = rdf.Filter("repetition == 0")
			.Filter("mwfrequence <= 1000")
			//.Filter("type != 2")
			.Histo1D({"Counts"," c to b",static_cast<int>(50),startPdf1, startPdf1 + 50*Params.FrequencyStep }, "mwfrequence");
	
	auto histF4 = 	rdf.Filter("repetition == 0")
			.Filter("mwfrequence >= 1000")
			//.Filter("type != 2")
			.Histo1D({"Counts"," d to a",static_cast<int>(50), startPdf2, startPdf2 + 50*Params.FrequencyStep}, "mwfrequence");

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
	histF4->SetTitle("d to a with cosmic events");
	histF4->SetMarkerStyle(21);
	histF4->SetMarkerColor(2);
	histF4->SetLineColor(4);
	histF4->DrawClone();
	histF4->DrawClone("SAME P");
        
	// Define some vectors to store the results
	vector<double> MCtruth;         // save MC true hyperfine splitting
	vector<double> v1_2017, v2_2017;// 
	vector<double> v1_rev, v2_rev;
	vector<double> v1_thr, v2_thr;
	vector<double> v1_cfrac, v2_cfrac;
	vector<double> v1_neigh, v2_neigh;
	vector<double> diff_2017, diff_rev, diff_thr, diff_cfrac, diff_neigh;
	
	//Onsetcb and Onset da
	vector<double> v1_cb, v2_da;
	
	int trial = 0;
	// LOOP ON TRIALS
	for(int i = 0; i < FileList.size(); i += 2){
		trial  += 1; std::cout << "Analizzo DataFrame " << trial << "\n" << std::endl;
		
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});
		
		// Extract from the dataset the starting frequency of the transition
		auto sweepStart_cb = frame.Filter("repetition == 0")
		                     .Filter("mwfrequence <= 1000")   // filter on the micro wave frequency to identify the c to b transition
		                     .Take<double>("sweepStart");
		auto actualStart_cb = sweepStart_cb.GetValue();       // get starting frequency from Rdataframe node
		
		// Load the data and save it in a Histogram (raw lineshape)
		// the histogram start from actualStart_cb, with a step given by Params.FrequencyStep
		auto SpectraCB = frame.Filter("repetition == 0")
				.Filter("mwfrequence <= 1000")   // filter on the micro wave frequency to identify the c to b transition
				 //.Filter("type != 2")          // uncomment to filter the background events
				.Histo1D({"Counts","mwfrequence", static_cast<int>(Params.SweepStep),actualStart_cb[0], actualStart_cb[0] + Params.SweepStep*Params.FrequencyStep }, "mwfrequence");
		
		
		// Extract from the dataset the starting frequency of the transition	
		auto sweepStart_da = frame.Filter("repetition == 0")
		                     .Filter("mwfrequence >= 1000")
		                     .Take<double>("sweepStart");
		auto actualStart_da = sweepStart_da.GetValue(); // get starting frequency from Rdataframe node
		
		// Load the data and save it in a Histogram (raw lineshape)
		 // the histogram start from actualStart_cb, with a step given by Params.FrequencyStep
		auto SpectraDA = frame.Filter("repetition == 0")
		                 .Filter("mwfrequence >= 1000")
		                 //.Filter("type != 2") // uncomment to filter the background events
		                 .Histo1D({"Counts","mwfrequence", static_cast<int>(Params.SweepStep), actualStart_da[0], actualStart_da[0] + Params.SweepStep*Params.FrequencyStep}, "mwfrequence");

		auto getOnsetCB = frame.Filter("repetition == 0").Filter("mwfrequence >= 1000").Take<double>("trueOnset");
		auto getOnsetDA = frame.Filter("repetition == 0").Filter("mwfrequence <= 1000").Take<double>("trueOnset");
		auto onsetda = getOnsetCB.GetValue(); //std::cout << "onset d to a : " << onsetda[0] << std::endl;
		auto onsetcb = getOnsetDA.GetValue(); //std::cout << "onset c to b : " << onsetcb[0] << std::endl;
		
		double onset1; // Reconstructed onset
		double onset2; // Reconstructed onset
		
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back((onsetda[0] - onsetcb[0]));
		
		// THRESHOLD
		onset1 = firstOverThreshold(SpectraCB, Nthr); // find onset cb
		onset2 = firstOverThreshold(SpectraDA, Nthr); // find onset da
		v1_thr.push_back(onset1 - (onsetcb[0]));
		v2_thr.push_back(onset2 - (onsetda[0]));
		diff_thr.push_back(onset2 - onset1 - (onsetda[0]  - onsetcb[0]));
		// 2017 FOWARD
		onset1 = algorithm_2017(SpectraCB); // find onset cb
		onset2 = algorithm_2017(SpectraDA); // find onset da
		
		v1_2017.push_back(onset1 - (onsetcb[0]));
		v2_2017.push_back(onset2 - (onsetda[0])); 
		diff_2017.push_back(onset2 - onset1 - (onsetda[0]   - onsetcb[0]));
		// REVERSED 2017 
		onset1 = reverse_2017(SpectraCB); // find onset cb
		onset2 = reverse_2017(SpectraDA); // find onset da
		
		v1_cb.push_back(onset1); v2_da.push_back(onset2);
		v1_rev.push_back(onset1 - (onsetcb[0]));
		v2_rev.push_back(onset2 - (onsetda[0]));
		diff_rev.push_back(onset2 - onset1 - (onsetda[0]   - onsetcb[0]));
		//std::cout << "onset 1 - onset 2 - MCtruth: " << onset2 - onset1 - (onsetda[0] -   - onsetcb[0]) << std::endl;
		// CONSTANT FRACTION
		onset1 = constFrac(SpectraCB, fraction, CosmicBackground,Nfilter); // find onset cb
		onset2 = constFrac(SpectraDA, fraction, CosmicBackground,Nfilter); // find onset da
		
		v1_cfrac.push_back(onset1 - (onsetcb[0]));
		v2_cfrac.push_back(onset2 - (onsetda[0]));
		diff_cfrac.push_back(onset2 - onset1 - (onsetda[0]  - onsetcb[0]));
		
		std::cout << std::setprecision(10);
		std::cout << "c to b onset: " <<  onset1 << std::endl;
		std::cout << "d to a onset: " <<  onset2 << std::endl;
		
		// SUM NEIGHBORS
		onset1 = Significance(SpectraCB, Nsigma, CosmicBackground, Nfilter); // find onset cb
		onset2 = Significance(SpectraDA, Nsigma, CosmicBackground, Nfilter); // find onset da
		
		v1_neigh.push_back(onset1 - (onsetcb[0]));
		v2_neigh.push_back(onset2 - (onsetda[0]));
		diff_neigh.push_back(onset2 - onset1 - (onsetda[0]  - onsetcb[0]));
		
	}
	
	// Plots
	auto h1 = new TH1D("h1","onset_{algorithm} - onset_{true}", 81, -40.5, 40.5);
	auto h2 = new TH1D("h2","onset_{algorithm} - onset_{true}", 81, -40.5 , 40.5);
	auto h3 = new TH1D("h3","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",61, -20.5 , 40.5);
	auto h4 = new TH1D("h4","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h5 = new TH1D("h5","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h6 = new TH1D("h6","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -300 , 300);
	auto h7 = new TH1D("h7","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h8 = new TH1D("h8","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h9 = new TH1D("h9","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	auto h10 = new TH1D("h10","onset_{algorithm} - onset_{true}", 121, -70.5, 70.5);
	auto h11 = new TH1D("h11","onset_{algorithm} - onset_{true}", 121, -70.5 , 70.5);
	auto h12 = new TH1D("h12","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 70.5);
	auto h13 = new TH1D("h13","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h14 = new TH1D("h14","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h15 = new TH1D("h15","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	//Plots for High rates
	auto h16 = new TH1D("h16","onset_{algorithm} c to b", 26, 147.5, 277.5);
	auto h17 = new TH1D("h17","onset_{algorithm} d to a", 26, 1420000 + 147.5 , 1420000 + 277.5);
	
	////////////////////////////////
	//2017 FOWARD
	std::vector<double> w(v1_2017.size(),1); // weights vector
	
   	h1->FillN(v1_2017.size(),v1_2017.data(), w.data());
   	h2->FillN(v2_2017.size(),v2_2017.data(), w.data());
   	h3->FillN(diff_2017.size(), diff_2017.data(), w.data());
   	
   	auto g = new TGraph(MCtruth.size(), MCtruth.data(), diff_2017.data());
	auto hh1 = new TH2D("hh1", "Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",100, 1419.96e3,1420.02e3,100,-11,+31); // 2d hist
   	hh1->FillN(MCtruth.size(), MCtruth.data(),diff_2017.data(),w.data());
   	
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
	pad2->cd(4);
	
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

 	h16->FillN(v1_cb.size(),v1_cb.data(), w.data());
   	h17->FillN(v2_da.size(),v2_da.data(), w.data());
	auto onsets = new TCanvas("onsets", "(onset_{cb} and  onset_{da})");
	auto pad7 = new TPad("pad7", "pad2",0,0,1,1);
	pad7->Divide(2,1,0.001,0.001);
	pad7->Draw();
	pad7->cd(1);
	h16->GetXaxis()->SetTitle("frequency [kHz]");
	h16->SetLineColor(2);
	h16->SetLineWidth(2);
	h16->Draw();
	pad7->cd(2);	
	h17->GetXaxis()->SetTitle("frequency [kHz]");
	h17->SetLineColor(2);
	h17->SetLineWidth(2);
	h17->Draw();
	
	////////////////////////////////
	// 2017 REVERSED
   	h4->FillN(v1_rev.size(),v1_rev.data(), w.data());
   	h5->FillN(v2_rev.size(),v2_rev.data(), w.data());
   	h6->FillN(diff_rev.size(), diff_rev.data(), w.data());
	
	auto hh2 = new TH2D("hh2", "Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",100, 1419.96e3,1420.02e3,100,-35,+35); // 2d hist
	hh2->FillN(MCtruth.size(), MCtruth.data(),diff_rev.data(),w.data());
	
	auto canvas1 = new TCanvas("d1", "2017 reversed", 1000,550);
	auto pad3 = new TPad("pad3", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad3->Divide(2,2,0.001,0.001);
	pad3->Draw();
	pad3->cd(1);
	h4->GetXaxis()->SetTitle("frequency [kHz]");
	h4->SetLineColor(2);
	h4->SetLineWidth(2);
	h4->Draw();
	pad3->cd(2);
	h5->GetXaxis()->SetTitle("frequency [kHz]");
	h5->SetLineColor(2);
	h5->SetLineWidth(2);
	h5->Draw();
	pad3->cd(3);
	h6->GetXaxis()->SetTitle("frequency [kHz]");
	h6->SetLineWidth(2);
	h6->SetLineColor(38);
	h6->Draw();
	name = TString::Format("2017_reversed");
	pad3->cd(4);
	hh2->Draw("COLZ");

	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas1->SaveAs(folder + name + endname);
	
	auto delta1 = new TCanvas("dd1", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h6->GetXaxis()->SetTitle("frequency [kHz]");
	h6->SetLineWidth(2);
	h6->SetLineColor(38);
	h6->Draw();
	name = "delta_2017_reversed";
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta1->SaveAs(folder + name + endname);
	
	////////////////////////////////
	// THRESHOLD
	h7->FillN(v1_thr.size(),v1_thr.data(), w.data());
   	h8->FillN(v2_thr.size(),v2_thr.data(), w.data());
   	h9->FillN(diff_thr.size(), diff_thr.data(), w.data());
	
	auto hh3 = new TH2D("hh3", "Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",100, 1419.96e3,1420.02e3,100,-35,+35); // 2d hist
	hh3->FillN(MCtruth.size(), MCtruth.data(),diff_thr.data(),w.data());
	
	auto canvas2 = new TCanvas("d2", "threshold", 1000,550);
	auto pad4 = new TPad("pad4", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad4->Divide(2,2,0.001,0.001);
	pad4->Draw();
	pad4->cd(1);
	h7->GetXaxis()->SetTitle("frequency [kHz]");
	h7->SetLineColor(2);
	h7->SetLineWidth(2);
	h7->Draw();
	pad4->cd(2);
	h8->GetXaxis()->SetTitle("frequency [kHz]");
	h8->SetLineColor(2);
	h8->SetLineWidth(2);
	h8->Draw();
	pad4->cd(3);
	h9->GetXaxis()->SetTitle("frequency [kHz]");
	h9->SetLineWidth(2);
	h9->SetLineColor(38);
	h9->Draw();
	pad4->cd(4);
	hh3->Draw("COLZ");
	
	name = TString::Format("Threshold_%d", static_cast<int>(Nthr));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas2->SaveAs(folder + name + endname);
	
	auto delta2 = new TCanvas("dd2", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h9->GetXaxis()->SetTitle("frequency [kHz]");
	h9->SetLineWidth(2);
	h9->SetLineColor(38);
	h9->Draw();
	name = TString::Format("delta_Threshold_%d", static_cast<int>(Nthr));
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta2->SaveAs(folder + name + endname);
	
	////////////////////////////////
	//COSTANT FRACTION
	h10->FillN(v1_cfrac.size(),v1_cfrac.data(), w.data());
   	h11->FillN(v2_cfrac.size(),v2_cfrac.data(), w.data());
   	h12->FillN(diff_cfrac.size(), diff_cfrac.data(), w.data());
   	auto g1 = new TGraph(MCtruth.size(), MCtruth.data(), diff_cfrac.data());
	
	auto hh4 = new TH2D("hh4", "Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",100, 1419.96e3,1420.02e3,100,-35,+35); // 2d hist
	hh4->FillN(MCtruth.size(), MCtruth.data(),diff_cfrac.data(),w.data());
	
	auto canvas3 = new TCanvas("d3", "constFract", 1000,550);
	auto pad5 = new TPad("pad5", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad5->Divide(2,2,0.001,0.001);
	pad5->Draw();
	pad5->cd(1);
	h10->GetXaxis()->SetTitle("frequency [kHz]");
	h10->SetLineColor(2);
	h10->SetLineWidth(2);
	h10->Draw();
	pad5->cd(2);
	h11->GetXaxis()->SetTitle("frequency [kHz]");
	h11->SetLineColor(2);
	h11->SetLineWidth(2);
	h11->Draw();
	pad5->cd(3);
	h12->GetXaxis()->SetTitle("frequency [kHz]");
	h12->SetLineWidth(2);
	h12->SetLineColor(38);
	h12->Draw();
	pad5->cd(4);
	hh4->Draw("COLZ");

	name = TString::Format("constFract_%d", static_cast<int>(100*fraction));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas3->SaveAs(folder + name + endname);
	
	auto delta3 = new TCanvas("dd3", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h12->GetXaxis()->SetTitle("frequency [kHz]");
	h12->SetLineWidth(2);
	h12->SetLineColor(38);
	h12->Draw();
	name = TString::Format("delta_constFract_%d", static_cast<int>(100*fraction));
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta3->SaveAs(folder + name + endname);

	////////////////////////////////
	// SIGNIFICANCE
	h13->FillN(v1_neigh.size(),v1_neigh.data(), w.data());
   	h14->FillN(v2_neigh.size(),v2_neigh.data(), w.data());
   	h15->FillN(diff_neigh.size(), diff_neigh.data(), w.data());
	auto g2 = new TGraph(MCtruth.size(), MCtruth.data(), diff_neigh.data());
	
	auto hh5 = new TH2D("hh5", "Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",100, 1419.96e3,1420.02e3,100,-35,+35); // 2d hist
	hh5->FillN(MCtruth.size(), MCtruth.data(),diff_neigh.data(),w.data());
	
	auto canvas4 = new TCanvas("d4", "sumNeighbours", 1000,550);
	auto pad6 = new TPad("pad6", "pad2",0,0,0,0);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad6->Divide(2,2,0.001,0.001);
	pad6->Draw();
	pad6->cd(1);
	h13->GetXaxis()->SetTitle("frequency [kHz]");
	h13->SetLineColor(2);
	h13->SetLineWidth(2);
	h13->Draw();
	pad6->cd(2);
	h14->GetXaxis()->SetTitle("frequency [kHz]");
	h14->SetLineColor(2);
	h14->SetLineWidth(2);
	h14->Draw();
	pad6->cd(3);
	h15->GetXaxis()->SetTitle("frequency [kHz]");
	h15->SetLineWidth(2);
	h15->SetLineColor(38);
	h15->Draw();
	pad6->cd(4);
	hh5->Draw("COLZ");
	/*
	g2->SetMarkerStyle(7); // Medium Dot
	g2->SetMarkerColor(9);
	g2->SetLineStyle(0);
	g2->SetTitle("Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth} ");
	g2->Draw("AP");
	*/
	name = TString::Format("significance_sigma%d", static_cast<int>(Nthr));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas4->SaveAs(folder + name + endname);
	
	auto delta4 = new TCanvas("dd4", "(onset_{pdf1} - onset_{pdf2}) - MC_{truth}");
	h15->GetXaxis()->SetTitle("frequency [kHz]");
	h15->SetLineWidth(2);
	h15->SetLineColor(38);
	h15->Draw();
	name = TString::Format("delta_significance_sigma%d", static_cast<int>(Nthr));
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	delta4->SaveAs(folder + name + endname);
	//TString name = TString::Format("sumNeighbours_sigma%d(cosmic=0)", static_cast<int>(Nthr)); TString endname = ".pdf"; 

}
