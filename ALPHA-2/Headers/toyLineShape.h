#ifndef FITALLFREQUENCY_H
#define FITALLFREQUENCY_H
#define VERBOSE 0
#include <iostream>
#include <string>
#include "TSystem.h"
#include "TH1.h"
#include "TSpline.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include <ROOT/RDataFrame.hxx>
#include "RooAbsReal.h"
#include "TString.h"
#include "TNtuple.h"
#include <TRandom3.h>

void SetContent(TH1 * histpdf, int Nbin, double (*f)(double, double, double, double), double start,double peak, double end){
	//USING A FUNCTION, EVALUATE THE FUNCTION AT X AND FILL HISTOGRAMS
	for(int i = 1; i <= Nbin; ++i){
		double x = histpdf->GetBinCenter(i);	// Point to be evaluated
		if(f(x, start,peak, end) > 0.){					// Check the function is > 0
			histpdf->SetBinContent(i,f(x, start, peak, end));
		} 
		else{
			histpdf->SetBinContent(i,0.);
		}// Set to 0
	}
	histpdf->Scale(1./histpdf->Integral(), "width");
}

void SetContent(TH1 * histpdf,
		int Nbin,
		double (*f)(double, double, double, double, double, double, double, double),
		double onset,
		double x0,
		double sigma0,
		double sigma1,
		double k0,
		double k1,
		double Norm){
	//USING A FUNCTION, EVALUATE THE FUNCTION AT X AND FILL HISTOGRAMS
	for(int i = 1; i <= Nbin; ++i){
		double x = histpdf->GetBinCenter(i);	// Point to be evaluated
		if(f(x,  onset, x0, sigma0, sigma1, k0, k1, Norm) > 0.){	// Check the function is > 0
			histpdf->SetBinContent(i,f(x, onset, x0, sigma0, sigma1, k0, k1, Norm));
		} 
		else{
			histpdf->SetBinContent(i,0.);
		}// Set to 0
	}
	histpdf->Scale(1./histpdf->Integral(), "width");
}

double ComputeProb(TH1 * histpdf, int i){
	// COMPUTE THE PROBABILITY OF THE BIN
	Double_t width = histpdf->GetBinWidth(i);
	Double_t prob = histpdf->GetBinContent(i);
	prob = prob*width;
	return prob;
}

void SetVectors(RooDataSet &Global,
		RooDataSet *data,
		vector<int> &b,
		int flag,
		int &Counts,
		vector<double> &f,
		double frequence,
		vector<int> &RunNumber,
		int run){
	// FILL VECTORS, DATASET...
	if(data){  
		Counts = data->sumEntries();
		Global.append(*data);
		for(int i = 0; i < data->sumEntries(); ++i){
		b.push_back(flag);
		f.push_back(frequence);
		RunNumber.push_back(run);
		}
	} else{ // If the dataset is Empy push back 0
	Counts = 0;
	}
}

void SetVectors(RooDataSet *data, vector<double> &a,vector<int> &b, int &Counts, vector<double> &f ,double frequence, int flag ){	// FILL VECTORS, DATASET...
	if(data){a.push_back(data->sumEntries());}
	Counts = data->sumEntries();
	for(int i = 0; i < data->sumEntries(); ++i){
	b.push_back(flag);
	f.push_back(frequence);
	}
}

void SetCoefficients(double a, double b, double c, RooRealVar *Nmix, RooRealVar *Ngas, RooRealVar *Nbk){
	Nmix->setVal(a);	// Set Nmix
	Ngas->setVal(b);	// Set Ngas
	Nbk->setVal(c); 
}

auto FillDataFrame(ROOT::RDataFrame &d1, RooDataSet &data,
				 	vector<double> &f,
				 	vector<int> &vType,
				 	vector<double> &vTot,
				 	int &j,
				 	vector<Double_t> &rn,
				 	vector<int> &RunNumber,
				 	vector<double> &freqDelay
				 	){
	
	auto rdf = d1
	.Define("runNumber", 
	[&RunNumber, &j](){
		return RunNumber[j];
	})
	.Define("random",
	[&rn, &j](){			//To subsample and randomize 
		return rn[j];
		})
	.Define("delay", [&freqDelay, &RunNumber, &j](){
		return freqDelay[RunNumber[j]];
	})
	.Define("frequence",	// Frequence of the event
	[&j, &f](){
		return f[j]; 
		})
	.Define("type", 		// Type of the event (0 wall, 1 res gas, 2 cosmic)
		[&j, &vType](){
		return vType[j];
		})
	.Define("radius",		// Generated radius
		[&data, &j, &vTot, &f](){
		if(VERBOSE){ 		// PRINT INFO ABOUT EVENTS
			std::cout << "Event id: " << j << " f: " << f[j] << std::endl; 
			data.get(j)->Print("V");
			}
		const RooArgSet &argSet = *(data.get(j));
		++j;				// Update event id
		return static_cast<RooAbsReal&>(argSet["x"]).getVal();
		});
	return rdf; // Return the node of the RDataFrame
}
#endif
