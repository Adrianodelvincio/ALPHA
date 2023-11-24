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

void ConvertTNtutpla(TNtuple &file_pdf, vector<double> &v1, vector<double> &v2){
	// CONVERT TNtuple IN A STD::VECTOR
	for (int row = 0; row < file_pdf.GetEntries(); ++row) {
	file_pdf.GetEntry(row);
	v1.push_back(file_pdf.GetArgs()[0]); 	// Extract frequence
	v2.push_back(file_pdf.GetArgs()[1]);	// Extract Nls
	}
}

void SetContent(TH1 * histpdf, int Nbin, double (*f)(double)){
	// IF USING SPLINE, EVALUATE THE SPLINE AT X AND FILL HISTOGRAMS
	for(int i = 1; i <= Nbin; ++i){
		if(f(histpdf->GetBinCenter(i)) > 0.){
			histpdf->SetBinContent(i,f(histpdf->GetBinCenter(i)));
		}
		else{
			histpdf->SetBinContent(i,0.);
		}
	}
}

void SetContent(TH1 * histpdf, int Nbin, TSpline3 * spline){
	// IF USING SPLINE, EVALUATE THE SPLINE AT X AND FILL HISTOGRAMS
	for(int i = 1; i <= Nbin; ++i){
		if(spline->Eval(histpdf->GetBinCenter(i)) > 0.){
		histpdf->SetBinContent(i,spline->Eval(histpdf->GetBinCenter(i)));}
		else{
		histpdf->SetBinContent(i,0.);}
		}
}

void SetNormalization(TH1 * histpdf){
	// NORMALIZE THE HISTOGRAM
	histpdf->Scale(1./histpdf->Integral(), "width");
}

double ComputeProb(TH1 * histpdf, int i){
	// COMPUTE THE PROBABILITY OF THE BIN
	Double_t width = histpdf->GetBinWidth(i);
	Double_t prob = histpdf->GetBinContent(i);
	prob = prob*width;
	return prob;
}

void SetVectors(RooDataSet &Global,RooDataSet *data, vector<double> &a,vector<int> &b, int &Counts, vector<double> &f ,double frequence, int flag ){
	// FILL VECTORS, DATASET...
	if(data){  
		a.push_back(data->sumEntries());
		Counts = data->sumEntries();
		Global.append(*data);
		}
	else{ // If the dataset is Empy push back 0
	a.push_back(0);
	Counts = 0;
	}
	for(int i = 0; i < data->sumEntries(); ++i){
	b.push_back(flag);
	f.push_back(frequence);
	}
}

void SetVectors(RooDataSet *data, vector<double> &a,vector<int> &b, int &Counts, vector<double> &f ,double frequence, int flag ){
	// FILL VECTORS, DATASET...
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

void FillDataFrame(ROOT::RDataFrame &d1, TString datafileName, TH1 * histpdf, RooDataSet data, vector<int> &vType, vector<double> &vTot, bool SaveData = true){
	
	int j(0); 	// Variable for loop
	int k(1); 	// Inner Loop, Events belonging to a single frequence
	int bin(1);	// Bin number 
	auto rdf = d1.Define("id", [&j](){		// Id of the events
		return j; 
		})
	.Define("frequence",	// Frequence of the event
	[&bin, &histpdf](){
		return histpdf->GetBinCenter(bin); 
		})
	.Define("Type", //
		[&vType, &j](){
		return vType[j];
		})
	.Define("radius",	// Generated radius
		[&j,&vTot, &data, &bin, &k](){
		if(k >=  vTot[bin-1]){
			++bin; 	// All counts per freq. are saved, update the bin
			k = 1;	// Set k to 0 for the next frequence inner loop
		}else{
			++k;	// Update inner loop
		}
		if(VERBOSE){
		std::cout << "Event id: " << j << std::endl; 
		std::cout << "Bin: " << bin-1 << " Counts per freq: " << vTot[bin -1] << " k event: " << k; 
		data.get(j)->Print("V");
		}
		const RooArgSet &argSet = *(data.get(j));
		++j;	// Update event id
		return static_cast<RooAbsReal&>(argSet["x"]).getVal();
		}); 
	rdf.Snapshot("myTree", datafileName);
}

auto FillDataFrame(ROOT::RDataFrame &d1, RooDataSet &data, vector<double> &f , vector<int> &vType, vector<double> &vTot,int &j){
	
	auto rdf = d1.Define("id",// Id of the events
	 [&j](){		
		return j;
		})
	.Define("frequence",// Frequence of the event
	[&j, &f](){
		return f[j]; 
		})
	.Define("Type", //
		[&j, &vType](){
		return vType[j];
		})
	.Define("radius",// Generated radius
		[&data, &j, &vTot, &f](){
		if(VERBOSE){ // PRINT INFO ABOUT EVENTS
			std::cout << "Event id: " << j << " f: " << f[j] << std::endl; 
			data.get(j)->Print("V");
			}
		const RooArgSet &argSet = *(data.get(j));
		++j;		// Update event id
		return static_cast<RooAbsReal&>(argSet["x"]).getVal();
		});
	return rdf;	// Return the node of the RDataFrame
}

#endif
