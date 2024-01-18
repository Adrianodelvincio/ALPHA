#include <iostream>
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "../Headers/AnalysisLineShape.h"
#include "../Headers/ConfigurationParser.h"
#include "../Headers/Algorithms.h"
using namespace RooFit;

double mean(std::vector<double> v){
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	return sum/ v.size();
}

double stdev(std::vector<double> v){
	double accum = 0.0;
	double m = mean(v);
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
    	accum += (d - m) * (d - m);
	});
	return sqrt(accum/(v.size() - 1));
}

std::vector<double> ScanAnalysis(TString directory,
					TString ConfFile,
					int start,
					int stop,
					double N,
					double fraction,
					double Nsigma,
					double Nfilter,
					double rate,
					TString folder
					){
	//ConfFile = directory + "ToyConfiguration.txt";
	//std::cout << ConfFile << std::endl;
	ReadConfFile Params(ConfFile); // Read the values from the configuration file 
	//Params.Print();
	Params.CosmicRate = rate;
	double CosmicBackground = Params.TimeStep * Params.CosmicRate;	// Number of Cosmic Events
	double FrequencyStep = Params.FrequencyStep;
	double startPdf1 = Params.x_cb_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep c-b
	double startPdf2 = Params.x_da_start - (FrequencyStep)*(Params.BinBeforeOnset + 0.5);	// Start of frequency sweep d-a
	double SweepStep = Params.SweepStep; // number of bin for each lineshape
	
	std::vector<std::string> FileList;
	FileList = getFiles(start,stop, directory); // file list to be analyzed

	// Define some vectors to store the results
	vector<double> MCtruth;
	vector<double> v1_2017, v2_2017;
	vector<double> v1bk_2017, v2bk_2017;
	vector<double> v1_rev, v2_rev;
	vector<double> v1bk_rev, v2bk_rev;
	vector<double> v1_thr, v2_thr;
	vector<double> v1bk_thr, v2bk_thr;
	vector<double> v1_cfrac, v2_cfrac;
	vector<double> v1_classic, v2_classic;
	vector<double> v1_neigh, v2_neigh;
	vector<double> v1_runningDiff, v2_runningDiff;
	vector<double> v1_hybrid, v2_hybrid;
	vector<double> diff_2017, diff_bk_2017;
	vector<double> diff_rev, diff_bk_rev;
	vector<double> diff_thr, diff_bk_thr;
	vector<double> diff_cfrac, diff_neigh;
	vector<double> cfrac_classic, diff_hybrid;
	vector<double> diff_runningDiff;
	int count = 0;
	// IMPLEMENTING THE TOY FOR THE ALGORITHM
	for(int i = 0; i < FileList.size(); i += 2){
		count  += 1; 
		std::cout << directory << " Analizzo DataFrame " << count << std::endl;
		
		ROOT::RDataFrame frame("myTree", {FileList[i], FileList[i+1]});		// Load i-th dataset
		
		/*
		auto Spectra1 = frame.Filter("frequence <= 1000")
					.Filter("type != 2 || runNumber == 1")
					.Histo1D({"Counts","Frequence", static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
		
		auto Spectra2 = frame.Filter("frequence >= 1000")
					.Filter("type != 2 || runNumber == 1")
					.Histo1D({"Counts","Frequence", static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");*/
		
		auto Spectra1 = frame.Filter("runNumber == 1")
			 .Filter("frequence <= 1000")
			 //.Filter("type != 2")
			 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep),startPdf1, startPdf1 + SweepStep*FrequencyStep }, "frequence");
		
		auto Spectra2 = frame.Filter("runNumber == 1")
			 .Filter("frequence >= 1000")
			 //.Filter("type != 2")
			 .Histo1D({"Counts","Frequence", static_cast<int>(SweepStep), startPdf2, startPdf2 + SweepStep*FrequencyStep}, "frequence");

		auto frame1 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Mean<double>("lineShift");
		auto frame2 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Mean<double>("lineShift");
		// Load the shifts of the lineshapes
		auto frame3 = frame.Filter("runNumber == 1").Filter("frequence >= 1000").Take<double>("lineShift");
		auto frame4 = frame.Filter("runNumber == 1").Filter("frequence <= 1000").Take<double>("lineShift");
		auto lineShiftda = frame3.GetValue();
		//std::cout << "LineShift d to a : " << lineShiftda[0] << std::endl;
		auto lineShiftcb = frame4.GetValue();
		//std::cout << "LineShift c to b : " << lineShiftcb[0] << std::endl;
		
		double onset1;				// Reconstructed onset
		double onset2;				// Reconstructed onset
		double threshold = N*CosmicBackground; // threshold considering the cosmic background
		
		// SAVE THE GENERATED MONTECARLO
		MCtruth.push_back((Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
		// THRESHOLD
			// With Cosmic Background
		onset1 = firstOverThreshold(Spectra1, N, CosmicBackground);
		onset2 = firstOverThreshold(Spectra2, N, CosmicBackground);
		v1bk_thr.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2bk_thr.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_bk_thr.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));		
			// Without Cosmic Background
		onset1 = firstOverThreshold(Spectra1, N);
		onset2 = firstOverThreshold(Spectra2, N);
		v1_thr.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_thr.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_thr.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
		// 2017 FOWARD
			// With Cosmic Background
		onset1 = algorithm_2017(Spectra1, CosmicBackground);
		onset2 = algorithm_2017(Spectra2, CosmicBackground);
		v1bk_2017.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2bk_2017.push_back(onset2 - (Params.x_da_start + lineShiftda[0])); 
		diff_bk_2017.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
			// Without Cosmic Background
		onset1 = algorithm_2017(Spectra1);
		onset2 = algorithm_2017(Spectra2);
		v1_2017.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_2017.push_back(onset2 - (Params.x_da_start + lineShiftda[0])); 
		diff_2017.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
		// REVERSED 2017
			// with cosmic background
		onset1 = reverse_2017(Spectra1, CosmicBackground);
		onset2 = reverse_2017(Spectra2, CosmicBackground);
		v1bk_rev.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2bk_rev.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_bk_rev.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
			// without cosmic background
		onset1 = reverse_2017(Spectra1);
		onset2 = reverse_2017(Spectra2);
		v1_rev.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_rev.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_rev.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
		// CONSTANT FRACTION
			// with cosmic subtraction
		onset1 = constFrac(Spectra1, fraction, CosmicBackground);
		onset2 = constFrac(Spectra2, fraction, CosmicBackground);
		v1_cfrac.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_cfrac.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_cfrac.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
			// without cosmic subtraction
		onset1 = constFrac(Spectra1, fraction);
		onset2 = constFrac(Spectra2, fraction);
		v1_classic.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_classic.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		cfrac_classic.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		
		// ALGORITMI CON FITRO/SOMMA SU BINS
		// SUM NEIGHBORS
		onset1 = sumNeighbors(Spectra1, Nsigma, CosmicBackground);
		onset2 = sumNeighbors(Spectra2, Nsigma, CosmicBackground);
		v1_neigh.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_neigh.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_neigh.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// Running Difference
		onset1 = runningDiff(Spectra1);
		onset2 = runningDiff(Spectra2);
		v1_runningDiff.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_runningDiff.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_runningDiff.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));
		// Hybrid Method
		onset1 = hybrid_cfSum(Spectra1, fraction ,CosmicBackground);
		onset2 = hybrid_cfSum(Spectra2, fraction ,CosmicBackground);
		v1_hybrid.push_back(onset1 - (Params.x_cb_start + lineShiftcb[0]));
		v2_hybrid.push_back(onset2 - (Params.x_da_start + lineShiftda[0]));
		diff_hybrid.push_back(onset2 - onset1 - (Params.x_da_start + lineShiftda[0] - Params.x_cb_start - lineShiftcb[0]));	
	}
	
	auto h1 = new TH1D("h1","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h2 = new TH1D("h2","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h3 = new TH1D("h3","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	auto h4 = new TH1D("h4","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h5 = new TH1D("h5","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h6 = new TH1D("h6","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	auto h7 = new TH1D("h7","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h8 = new TH1D("h8","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h9 = new TH1D("h9","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	auto h10 = new TH1D("h10","onset_{algorithm} - onset_{true}", 121, -70.5, 50.5);
	auto h11 = new TH1D("h11","onset_{algorithm} - onset_{true}", 121, -70.5 , 50.5);
	auto h12 = new TH1D("h12","(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",121, -70.5 , 50.5);
	
	TString name;
	TString endname = ".pdf"; 
	int numero = 0; 
	std::vector<double> w(v1_2017.size(),1); // weights vector
	
	//2017 FOWARD
	
   	h10->FillN(v1_2017.size(),v1_2017.data(), w.data());
   	h11->FillN(v2_2017.size(),v2_2017.data(), w.data());
   	h12->FillN(diff_2017.size(), diff_2017.data(), w.data());
   	
   	auto g = new TGraph(MCtruth.size(), MCtruth.data(), diff_2017.data());
	auto hh10 = new TH2D("hh10", "Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth}",100, 1419.96e3,1420.02e3,100,-35,+35); // 2d hist
   	hh10->FillN(MCtruth.size(), MCtruth.data(),diff_2017.data(),w.data());
   	
   	auto canvas = new TCanvas("d", "2017 algorithm", 1000,550);
	auto pad2 = new TPad("pad2", "pad2",0,0,1,1);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad2->Divide(2,2,0.001,0.001);
	pad2->Draw();
	pad2->cd(1);
	h10->GetXaxis()->SetTitle("frequency [kHz]");
	h10->SetLineColor(2);
	h10->SetLineWidth(2);
	h10->Draw();
	pad2->cd(2);
	h11->GetXaxis()->SetTitle("frequency [kHz]");
	h11->SetLineColor(2);
	h11->SetLineWidth(2);
	h11->Draw();
	pad2->cd(3);
	h12->GetXaxis()->SetTitle("frequency [kHz]");
	h12->SetLineWidth(2);
	h12->SetLineColor(38);
	h12->Draw();
	pad2->cd(4);
	hh10->Draw("COLZ");
	
	name = "2017_foward";
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas->SaveAs(folder + name + endname);
	
	//COSTANT FRACTION
	h1->FillN(v1_cfrac.size(),v1_cfrac.data(), w.data());
   	h2->FillN(v2_cfrac.size(),v2_cfrac.data(), w.data());
   	h3->FillN(diff_cfrac.size(), diff_cfrac.data(), w.data());
   	auto g1 = new TGraph(MCtruth.size(), MCtruth.data(), diff_cfrac.data());
	
	auto canvas3 = new TCanvas("d3", "constFract", 1000,550);
	auto pad5 = new TPad("pad5", "pad2",0,0,1,1);
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
	pad5->cd(4);
	g1->SetMarkerStyle(7); // Medium Dot
	g1->SetMarkerColor(9);
	g1->SetLineStyle(0);
	g1->SetTitle("Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth} ");
	g1->Draw("AP");
	
	name = TString::Format("constFract_%d", static_cast<int>(100*fraction));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas3->SaveAs(folder + name + endname);

	// SUMNEIGHBORS	
	h4->FillN(v1_neigh.size(),v1_neigh.data(), w.data());
   	h5->FillN(v2_neigh.size(),v2_neigh.data(), w.data());
   	h6->FillN(diff_neigh.size(), diff_neigh.data(), w.data());
	auto g2 = new TGraph(MCtruth.size(), MCtruth.data(), diff_neigh.data());
	
	auto canvas4 = new TCanvas("d4", "sumNeighbours", 1000,550);
	auto pad6 = new TPad("pad6", "pad2",0,0,1,1);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad6->Divide(2,2,0.001,0.001);
	pad6->Draw();
	pad6->cd(1);
	h4->GetXaxis()->SetTitle("frequency [kHz]");
	h4->SetLineColor(2);
	h4->SetLineWidth(2);
	h4->Draw();
	pad6->cd(2);
	h5->GetXaxis()->SetTitle("frequency [kHz]");
	h5->SetLineColor(2);
	h5->SetLineWidth(2);
	h5->Draw();
	pad6->cd(3);
	h6->GetXaxis()->SetTitle("frequency [kHz]");
	h6->SetLineWidth(2);
	h6->SetLineColor(38);
	h6->Draw();
	pad6->cd(4);
	g2->SetMarkerStyle(7); // Medium Dot
	g2->SetMarkerColor(9);
	g2->SetLineStyle(0);
	g2->SetTitle("Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth} ");
	g2->Draw("AP");
	
	name = TString::Format("sumNeighbors_sigma%d", static_cast<int>(Nsigma));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas4->SaveAs(folder + name + endname);
	
	//Hybrid
	h7->FillN(v1_hybrid.size(),v1_hybrid.data(), w.data());
   	h8->FillN(v2_hybrid.size(),v2_hybrid.data(), w.data());
   	h9->FillN(diff_hybrid.size(), diff_hybrid.data(), w.data());
   	auto g3 = new TGraph(MCtruth.size(), MCtruth.data(), diff_hybrid.data());
	
	auto canvas5 = new TCanvas("d5", "Hybrid", 1000,550);
	auto pad7 = new TPad("pad7", "pad2",0,0,1,1);
	gStyle->SetOptStat(1);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	pad7->Divide(2,2,0.001,0.001);
	pad7->Draw();
	pad7->cd(1);
	h7->GetXaxis()->SetTitle("frequency [kHz]");
	h7->SetLineColor(2);
	h7->SetLineWidth(2);
	h7->Draw();
	pad7->cd(2);
	h8->GetXaxis()->SetTitle("frequency [kHz]");
	h8->SetLineColor(2);
	h8->SetLineWidth(2);
	h8->Draw();
	pad7->cd(3);
	h9->GetXaxis()->SetTitle("frequency [kHz]");
	h9->SetLineWidth(2);
	h9->SetLineColor(38);
	h9->Draw();
	pad7->cd(4);
	g3->SetMarkerStyle(7); // Medium Dot
	g3->SetMarkerColor(9);
	g3->SetLineStyle(0);
	g3->SetTitle("Scatter plot; MC_{truth} ;(onset_{pdf1} - onset_{pdf2}) - MC_{truth} ");
	g3->Draw("AP");
	
	name = TString::Format("hybrid_%d", static_cast<int>(100*fraction));
	numero = 0;
	while(!gSystem->AccessPathName(folder + name + endname)){
		numero += 1;
		TString add = TString::Format("_%d", numero);
		name = name + add;
	}
	canvas5->SaveAs(folder + name + endname);
	
	//delete h1,h2,h3,h4,h5,h6;
	
	//std::cout << "rate " << rate << " reversed variance: " << stdev(diff_rev) << " foward variance: " << stdev(diff_2017) << std::endl;
	return 		{mean(diff_thr),	stdev(diff_thr),         // THRESHOLD
			mean(diff_bk_thr),	stdev(diff_bk_thr),		//with background
			mean(diff_2017),	stdev(diff_2017),        // FOWARD
			mean(diff_bk_2017),	stdev(diff_bk_2017),		//with background
			mean(diff_rev),		stdev(diff_rev),         // REVERSED
			mean(diff_bk_rev),	stdev(diff_bk_rev),		//with background
			mean(cfrac_classic),	stdev(cfrac_classic), 	// CONSTANT FRACTION 
			mean(diff_cfrac),	stdev(diff_cfrac),       	//with backgrounf
			mean(diff_neigh),	stdev(diff_neigh),       // SUM NEIGHBORS
			mean(diff_hybrid),	stdev(diff_hybrid),      // HYBRID METHOD
			mean(v1_thr),		mean(v1bk_thr),
			mean(v1_2017),		mean(v1bk_2017),
			mean(v1_rev), 		mean(v1bk_rev),
			mean(v1_classic),	mean(v1_cfrac),
			mean(v1_neigh),		mean(v1_hybrid),
			mean(v2_thr),		mean(v2bk_thr),
			mean(v2_2017),		mean(v2bk_2017),
			mean(v2_rev), 		mean(v2bk_rev),
			mean(v2_classic),	mean(v2_cfrac),
			mean(v2_neigh),		mean(v2_hybrid),
			mean(diff_runningDiff),	stdev(diff_runningDiff)};
}

