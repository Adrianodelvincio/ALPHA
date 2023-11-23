#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "Headers/toyLineShape.h"

using namespace RooFit;

void LoopLineShape(double Mix_c = 0.5, double Mix_d = 0.5, double C = 0.5, int NBin = 30, int NTOT = 10000, int Nloop = 1, bool Save = false){
	/* Parameters of the Simulation */
int Nbin = NBin;		// Number of Bins
int Ntot = NTOT;		// Number of Total Events
double Ncosmic = (0.492 * Nbin);// Number of Cosmic Events
double pWall_c = Mix_c;	// Weight annihilation on walls for pdf1 (transition c -> b)
double pWall_d = Mix_d;	// Weight annihilation on walls for pdf2 (transition d -> a)
double c = C;			// Percentage of division two datasets
	/*				*/
int Ntrial = Nloop;
double d = 1 - c; double Nc = Ntot*c; double Nd = Ntot*d;
double pGas_d = 1 - pWall_d; double pGas_c = 1 - pWall_c;

gInterpreter->GenerateDictionary("ToyLine", "Headers/toyLineShape.h");
TNtuple file_pdf1("pdf1", "pdf1","x:y");
TNtuple file_pdf2("pdf2", "pdf2","x:y");
file_pdf1.ReadFile("LineShape/lineShape1.csv");
file_pdf2.ReadFile("LineShape/lineShape2.csv");

vector<double> v1; // Frequency pd1
vector<double> v2; // Counts per frequence pdf1
vector<double> t1; // Frequence pdf2
vector<double> t2; // Counts per frequence pdf2

ConvertTNtutpla(file_pdf1,v1,v2);
ConvertTNtutpla(file_pdf2,t1,t2);

Double_t frequence[v1.size()]; Double_t pdf1[v2.size()];
Double_t frequence2[t1.size()]; Double_t pdf2[t2.size()];

std::copy(v1.begin(),v1.end(),frequence); std::copy(v2.begin(),v2.end(),pdf1);
std::copy(t1.begin(),t1.end(),frequence2); std::copy(t2.begin(),t2.end(),pdf2);

// Interpolate the data with spline
TSpline3 *spline1 = new TSpline3("LineShape1", frequence, pdf1, v1.size());
TSpline3 *spline2 = new TSpline3("LineShape2", frequence2, pdf2, t1.size());
TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin,frequence[0], frequence[v1.size() -1]);
TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin,frequence2[0], frequence2[t1.size() -1]);

// Set Content Histograms
SetContent(histpdf1,Nbin,spline1);
SetContent(histpdf2,Nbin,spline2);
// Normalize histograms
Double_t factor = 1.;
histpdf1->Scale(factor/histpdf1->GetEntries());
histpdf2->Scale(factor/histpdf2->GetEntries());

RooRealVar x("x","r [cm]",0.,4.);
//x.setBins(Nbin);

RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);

//MIXING analytic model
RooRealVar mu("mu", "mu", 2.38);
RooRealVar sigMix("sigMix", "sigMix",0.85);
RooGaussian gauss_Mix("gauss", "gauss", x, mu, sigMix);

//PDF UWlosses
//FIXING sigRay to the value from the fit to the model (see analysisMLE.cpp)
RooRealVar sigRay("sigRay", "sigma", 1.722);
RooGenericPdf Rayleigh("line", "linear model", " TMath::Abs(x)/(sigRay*sigRay) * TMath::Exp(-(x*x)/(2*sigRay*sigRay))", RooArgSet(x,sigRay));

//PDF Cosmic
//RooRealVar q("q", "q",0.25, 0.25 - 200);
RooGenericPdf linearFit("linearFit", "linear model", "TMath::Abs((0.125)*x)", RooArgSet(x));

RooRealVar Nmix("Nmix","Nmix",100, -3000, +1000000);
RooRealVar Ngas ("Ngas", "Ngas",100, -3000, +3000);
RooRealVar Nbk ("Ncosmic", "Ncosmic", 100, -3000, +3000);
//Model to generate the data
RooAddPdf genMix("model", "model", RooArgList{gauss_Mix}, RooArgList{Nmix});
RooAddPdf genGas("model1", "model1", RooArgList{Rayleigh}, RooArgList{Ngas});
RooAddPdf genCosmic("model3", "model3", RooArgList{linearFit}, RooArgList{Nbk});

vector<double> f1; 	// Frequencies pdf1
vector<double> v1Nmix;	// Counts Nmix 
vector<double> v1Ngas;	// Counts Ngas
vector<double> v1Nbk;	// Cosmic Events
vector<double> v1Tot;	// Total Counts
vector<int>    v1Type;	// Type of Event

vector<double> f2;	// Frequencies pdf2
vector<double> v2Nmix;  // Counts Nmix
vector<double> v2Ngas;	// Counts Ngas
vector<double> v2Nbk;	// Cosmic events
vector<double> v2Tot;	// Total Counts
vector<int>    v2Type;	// Type of Event

RooDataSet data("data", "data", RooArgSet(x)); // Dataset to store the events
RooDataSet dati("dati", "dati", RooArgSet(x));

//External Toy Loop
	for(int j = 0; j < Ntrial; j++){
		
		// LOOP, Assign the Nmix and Ngas per frequence
		int trueTot1 = 0; int trueTot2 = 0;
		for(int i = 1; i <= Nbin; i++){
			//Pdf1
			double prob = ComputeProb(histpdf1,i);	// Probability of the bin
			SetCoefficients((pWall_c*Nc)*prob,(pGas_c*Nc)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			int gasCount = 0; int mixCount = 0; int CosmicCount = 0;
			// MIX
			RooDataSet *dataLoopWall = genMix.generate(x,Extended());
			if(dataLoopWall){			// If Nmix different from 0
				data.append(*dataLoopWall);	// Append to big dataset
				SetVectors(dataLoopWall, v1Nmix, v1Type, mixCount, 0);
			}else {v1Nmix.push_back(0);}
			//GAS 
			RooDataSet *dataLoopGas = genGas.generate(x,Extended());
			if(dataLoopGas){			// If Ngas different from 0
				data.append(*dataLoopGas);
				SetVectors(dataLoopGas,v1Ngas,v1Type,gasCount,1);
			}else{ v1Ngas.push_back(0);}
			//COSMIC
			RooDataSet *dataCosmic = genCosmic.generate(x, Extended());
			if(dataCosmic){
				data.append(*dataCosmic);
				SetVectors(dataCosmic,v1Nbk,v1Type,CosmicCount,2);
			}else{v1Nbk.push_back(0);}
			
			
			f1.push_back(histpdf1->GetBinCenter(i)); 	// save frequency
			v1Tot.push_back(mixCount + gasCount + CosmicCount);
			trueTot1 += mixCount + gasCount + CosmicCount;
			
			//Pdf2
			prob = ComputeProb(histpdf2,i);
			SetCoefficients((pWall_d*Nd)*prob,(pGas_d*Nd)/Nbin,Ncosmic/Nbin, &Nmix,&Ngas,&Nbk);
			// MIX
			RooDataSet *dataLoopWall2 = genMix.generate(x,Extended());// Generate data
			if(dataLoopWall){	// Append to big dataset
				dati.append(*dataLoopWall2);		
				SetVectors(dataLoopWall2, v2Nmix, v2Type, mixCount, 0);
			}else{ v2Nmix.push_back(0);}
			//GAS
			RooDataSet *dataLoopGas2 = genGas.generate(x,Extended());
			if(dataLoopGas){	// If Ngas different from 0
				dati.append(*dataLoopGas2);
				SetVectors(dataLoopGas2,v2Ngas,v2Type,gasCount,1);
			}else{ v2Ngas.push_back(0);}
			//COSMIC
			if(dataCosmic){
				dati.append(*dataCosmic);
				SetVectors(dataCosmic,v2Nbk,v2Type,CosmicCount,2);
			}else{v2Nbk.push_back(0);}
			f2.push_back(histpdf2->GetBinCenter(i));
			v2Tot.push_back(mixCount + gasCount + CosmicCount);
			trueTot2 += mixCount + gasCount + CosmicCount;
		
		}

	ROOT::RDataFrame d1(trueTot-1); // PDF1
	ROOT::RDataFrame d2(trueTot -1); // PDF2

	// 
	int j(0); 	// Variable for loop
	int k(1); 	// Inner Loop, Events belonging to a single frequence
	int bin(1);	// Bin number 
	TString datafileName = TString::Format("LineShape/ToyShape1_%d_%d_c%d.root", static_cast<int>(pWall_c*100),
	static_cast<int>(pWall_d*100),
	static_cast<int>(c*100));

	d1.Define("id", [&j](){return j;})
		.Define("frequence",[&bin, &histpdf1](){return histpdf1->GetBinCenter(bin);})
		.Define("Type",[&v1Type, &j](){return v1Type[j];})
		.Define("radius",[&j,&v1Tot, &data, &bin, &k](){
			if(k >=  v1Tot[bin-1]){
				++bin; 	// All counts per freq. are saved, update the bin
				k = 1;	// Set k to 0 for the next frequence inner loop
			}else{
				++k;	// Update inner loop
			}
			const RooArgSet &argSet = *(data.get(j));
			++j;		// Update event id
			return static_cast<RooAbsReal&>(argSet["x"]).getVal();
			});

	TString datafile = TString::Format("LineShape/ToyShape2_%d_%d_c%d.root", static_cast<int>(pWall_c*100),
	static_cast<int>(pWall_d*100),
	static_cast<int>(c*100));

	j = 0; k = 1; bin = 1;
	d2.Define("id", [&j](){return j;})
		.Define("frequence", [&bin, &histpdf2](){ return histpdf2->GetBinCenter(bin);})
		.Define("Type", [&v2Type, &j](){return v2Type[j];})
		.Define("radius",[&j,&v2Tot, &dati, &bin, &k](){
			if(k >=  v2Tot[bin-1]){
				++bin; 	// All counts per freq. are saved, update the bin
				k = 1;	// Set k to 0 for the next frequence inner loop
			}else{
				++k;	// Update inner loop
			}
			const RooArgSet &argSet = *(dati.get(j));
			++j;		// Update event id
			return static_cast<RooAbsReal&>(argSet["x"]).getVal();
			});
		
		} // Loop Trials
	
} // End program

