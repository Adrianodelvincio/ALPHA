#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "../Headers/toyLineShape.h"

using namespace RooFit;

double modello1(double x, double xmin = 0){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= 1){ return 2*x;}
	else {return 0;}
}

double modello2(double x, double xmin = 1420){
	double xmax = 1422;
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return 2*(x - xmin);}
	else {return 0;}
}

double parabola1(double x, double xmin = 0){
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= 1){ return 3*(x*x);}
	else {return 0;}
}

double parabola2(double x, double xmin = 1420){
	double xmax = 1422;
	if(x <= xmin){ return 0;}
	else if( x > xmin && x <= xmax){ return (3/2 )*(x - xmin)*(x - xmin);}
	else {return 0;}
}

double algorithm(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);
double doubleThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);

void LoopLineShape(int Nloop = 400, double Mix_c = 1, double Mix_d = 1, double C = 0.5, int NBin = 30, int NTOT = 400, double Rising1 = 0, double Rising2 = 1420 ,bool Save = false, bool MethodSpline = false){
	/* Parameters of the Simulation */
int Nbin = NBin;		// Number of Bins
int Ntot = NTOT;		// Number of Total Events
double pWall_c = Mix_c;	// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Mix_d;	// Weight annihilation on walls for pdf2 (transition d -> a)
double c = C;			// Percentage of division two datasets
	/*				*/

int Ntrial = Nloop;
double d = 1 - c; double Nc = Ntot*c; double Nd = Ntot*d;
double pGas_d = 1 - pWall_d; double pGas_c = 1 - pWall_c;
double timeStep = 8; // time duration for each step
double Ncosmic = Nbin * timeStep * (1786./35000.);// Number of Cosmic Events

gInterpreter->GenerateDictionary("ToyLine","../Headers/toyLineShape.h");

//Histogram to store the distributions
double rangepdf1[2] = {-1.2096, 1.80480};
double rangepdf2[2] = {1419.20962, 1423.86533};
TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin, rangepdf1[0], rangepdf1[1]);
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin, rangepdf2[0], rangepdf2[1] );

if(MethodSpline){
	SplineMethod(histpdf1,histpdf2, Nbin);
}

if(!MethodSpline){
	SetContent(histpdf1,Nbin,parabola1, Rising1);
	SetContent(histpdf2,Nbin,parabola2, Rising2);
}

SetNormalization(histpdf1);
SetNormalization(histpdf2);

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

RooRealVar x("x","r [cm]",0.,4.);
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix); //Pdf Annihilation on walls
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay)); //PDF Residual Gas
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x)); //PDF Cosmic

RooRealVar Nmix("Nmix","Nmix",100, -3000, +1000000);
RooRealVar Ngas ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk ("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model", "model", RooArgList{gauss_Mix}, RooArgList{Nmix});
RooAddPdf genGas("model1", "model1", RooArgList{Rayleigh}, RooArgList{Ngas});
RooAddPdf genCosmic("model2", "model2", RooArgList{linearFit}, RooArgList{Nbk});

vector<double> onset1v, onset2v;  // Reconstructed onsets
vector<double> deltaOnset;

//External Toy Loop
	for(int l = 0; l < Ntrial; l++){
		RooDataSet dataPdf1("data", "data", RooArgSet(x)); // Dataset to store the events
		RooDataSet dataPdf2("dati", "dati", RooArgSet(x));
		vector<double> f1, f2;			// Frequencies pdf1 , Frequencies pdf2
		vector<double> v1Tot, v2Tot;	// Total Counts
		vector<int>    v1Type,v2Type;	// Type of Event
		
		for(int i = 1; i <= Nbin; i++){ // Inner loop, assign counts to each frequency
			//Pdf1
			double prob = ComputeProb(histpdf1,i);	// Probability of the bin
			SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
			
			if(prob > 0){
			RooDataSet *dataLoopWall = genMix.generate(x,Extended()); // GENERATE THE DATA
			SetVectors(dataPdf1,dataLoopWall, v1Type, mixCount,  f1,histpdf1->GetBinCenter(i),0); // FILL DATASET
			}else { mixCount = 0;}
			
			if(Ngas.getVal() != 0){
			RooDataSet *dataLoopGas = genGas.generate(x, Extended());
			SetVectors(dataPdf1,dataLoopGas, v1Type, gasCount,   f1,histpdf1->GetBinCenter(i),1);
			}else {gasCount = 0;}
			
			RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
			SetVectors(dataPdf1,dataCosmic,  v1Type, CosmicCount,f1,histpdf1->GetBinCenter(i),2);
			
			// STORE SOME USEFUL QUANTITIES
			//std::cout << "conteggi residual gas: " << gasCount << std::endl;
			//std::cout << "PROB COSMIC " << static_cast<double>(Ncosmic/Nbin) << std::endl;
			v1Tot.push_back(mixCount + gasCount + CosmicCount);
			
			// PDF 2
			mixCount = 0; gasCount = 0; CosmicCount = 0;
			prob = ComputeProb(histpdf2,i);	// Probability of the bin
			SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			
			if(prob > 0){
			RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());	// Generate Wall data
			SetVectors(dataPdf2,dataLoopWall2, v2Type, mixCount, f2,histpdf2->GetBinCenter(i),0);
			} else{  mixCount = 0;}
			if(Ngas.getVal() != 0){
			RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());	// Generate Gas counts
			SetVectors(dataPdf2,dataLoopGas2, v2Type, gasCount,   f2,histpdf2->GetBinCenter(i),1);
			} else{ gasCount = 0;}
			RooDataSet *dataCosmic2 = genCosmic.generate(x, Extended());			
			SetVectors(dataPdf2,dataCosmic2,   v2Type, CosmicCount,f2,histpdf2->GetBinCenter(i),2);
			
			// STORE USEFUL QUANTITIES
			//std::cout << "conteggi residual gas: " << gasCount << std::endl;
			v2Tot.push_back(mixCount + gasCount + CosmicCount);
		}
	
		int Tot1 = std::accumulate(v1Tot.begin(), v1Tot.end(), 0);
		int Tot2 = std::accumulate(v2Tot.begin(), v2Tot.end(), 0);
		
		// DATA ANALYSIS ON SIMULATED EVENTS
		//ROOT::RDataFrame d1(Tot1 -1); // PDF1
		ROOT::RDataFrame d1(Tot1 -1);
		ROOT::RDataFrame d2(Tot2 -1); // PDF2
		
		int j(0); // FROM HERE, APPLY THE ALGORITHM TO PDF1
		auto FilledFrame1 = FillDataFrame(d1, dataPdf1, f1, v1Type, v1Tot,j); // Fill RDataFrame
		auto histF1 = FilledFrame1.Histo1D({"Counts","Frequence", 60u,rangepdf1[0], rangepdf1[1]}, "frequence");
		auto ww = FilledFrame1.Display({"frequence", "Type"}, 1); ww->Print();
		
		double onset1; int bin;			// Reconstructed onset and bin of the onset
		double threshold = 5*Ncosmic/Nbin;// threshold considering the cosmic background
		onset1 = doubleThreshold(histF1, threshold); 
		onset1v.push_back(onset1 - Rising1);
		
		if(l == 0){ //save dataframe
			j= 0;
			std::cout << "Save LoopData1.root" << std::endl;
			auto rdf = FilledFrame1;
			rdf.Snapshot("myTree", "LoopData_1.root");
		}

		j = 0; // FROM HERE APPLY THE ALGORITHM TO PDF2
		auto FilledFrame2 = FillDataFrame(d2, dataPdf2, f2, v2Type, v2Tot,j); // Fill RDataFrame
		auto histF2 = FilledFrame2.Histo1D({"Counts","Frequence", 60u,rangepdf2[0], rangepdf2[1]}, "frequence");
		auto hh = FilledFrame2.Display({"frequence", "Type", "radius"}, 1); hh->Print();
		
		double onset2; bin = 0;			// Reconstructed onset and bin of the onset
		threshold = 5*Ncosmic/Nbin;		// threshold considering the cosmic background
		onset2 = doubleThreshold(histF2, threshold); 
		onset2v.push_back(onset2 - Rising2);
		
		deltaOnset.push_back((onset2 - onset1) - (Rising2 - Rising1));
		
		if(l == 0){ // Save dataframe
			j = 0;
			std::cout << "Save LoopData2.root" << std::endl;
			auto rdf = FilledFrame2;
			rdf.Snapshot("myTree", "LoopData_2.root");
		}
		//std::cout << "Total Event Gen. pdf1: " << Tot1 << std::endl;
		//std::cout << "Total Event Gen. pdf2: " << Tot2 << std::endl;
		std::cout << "Find onset with threshold: " << threshold <<std::endl;
		std::cout << "	Onset pdf1 : " << onset1 << std::endl;
		std::cout << "	Onset pdf2 : " << onset2 << std::endl;
		std::cout << "	delta Onset - MC truth : " << deltaOnset[l] << std::endl; 
	} // Loop Trials
	
	// Create the histograms
	std::vector<double> w(onset1v.size(),1); // weights vector
	auto h1 = new TH1D("h1","onset_{algorithm} - onset_{true}",20,-1.5,1);
   	h1->FillN(onset1v.size(),onset1v.data(), w.data());
   	auto h2 = new TH1D("h2","onset_{algorithm} - onset_{true}",20,-1.5,1);
   	h2->FillN(onset2v.size(),onset2v.data(), w.data());
   	auto h3 = new TH1D("h3","(onset_{1} - onset_{2}) - MC truth",20,-1.5,1);
   	h3->FillN(deltaOnset.size(),deltaOnset.data(), w.data());
   	
   	auto canvas = new TCanvas("d", "Toy Result", 1000,550);
	auto pad = new TPad("pad", "pad",0,0,1,1);
	gStyle->SetOptStat(1);
	//gStyle->SetOptFit(1);
	pad->Divide(2,1,0.005,0.005);
	pad->Draw();
	pad->cd(1);
	h1->SetLineColor(1);
	h1->SetLineWidth(2);
	h1->Draw();
	//legend->Draw();
	pad->cd(2);
	h2->SetLineColor(1);
	h2->SetLineWidth(2);
	h2->Draw();
	//legend->Draw();
	canvas->SaveAs("ToyResult_doubleThreshold.pdf");
	 	
	auto canvas1 = new TCanvas("d1", "Toy Result", 1000,550);
	h3->SetLineColor(3);
	h3->Draw(); 
} // End program

double algorithm(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) >= threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
	}
	std::cout << "		i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double doubleThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
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
