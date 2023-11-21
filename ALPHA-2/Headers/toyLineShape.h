#ifndef FITALLFREQUENCY_H
#define FITALLFREQUENCY_H
#include <iostream>
#include <string>
#include "TSystem.h"
#include "TH1.h"
#include "TSpline.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

void SetContent(TH1 * histpdf, int Nbin, TSpline3 * spline){
	for(int i = 1; i <= Nbin; ++i){
		if(spline->Eval(histpdf->GetBinCenter(i)) > 0.){
		histpdf->SetBinContent(i,spline->Eval(histpdf->GetBinCenter(i)));}
		else{
		histpdf->SetBinContent(i,0.);}
		}
}

double ComputeProb(TH1 * histpdf, int i){
	Double_t width = histpdf->GetBinWidth(i);
	Double_t prob = histpdf->GetBinContent(i);
	prob = prob*width;
	return prob;
}

void SetVectors(RooDataSet *data, vector<double> &a,vector<int> &b, int &Counts , int flag ){
	if(data){a.push_back(data->sumEntries());}
	Counts = data->sumEntries();
	for(int i = 0; i < data->sumEntries(); ++i){
	b.push_back(flag);
	}
}

void SetCoefficients(double a, double b, double c, RooRealVar *Nmix, RooRealVar *Ngas, RooRealVar *Nbk){
	Nmix->setVal(a);	// Set Nmix
	Ngas->setVal(b);	// Set Ngas
	Nbk->setVal(c); 
}
#endif
