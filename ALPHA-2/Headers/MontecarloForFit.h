#ifndef MONTECARLOFORFIT_H
#define MONTECARLOFORFIT_H
#include <iostream>
#include "RooRealVar.h"
#include <numeric>

void SetVectors(int i,std::vector<double> &mix,std::vector<double> &uw, std::vector<double> &bk, RooRealVar a,RooRealVar b,RooRealVar c,RooRealVar aa,RooRealVar bb,RooRealVar cc){
		mix.push_back(a.getVal() - aa.getVal());
		uw.push_back(b.getVal() - bb.getVal());
		bk.push_back(c.getVal() - cc.getVal());}
		

Double_t media(std::vector<Double_t>&v) {
if(v.size()<=0.) return 0.;
return (Double_t) std::accumulate(v.begin(), v.end(), 0.0) / (Double_t)v.size();
}

void SetProb(double *a, double* b, double* c, int N, int Nfix){
*c = Nfix/N;	// c fixed
double tt = *a;
*a = *a/(*a +*b); *b = *b/(tt+*b);
*a = (N - Nfix)*(*a/N); // Correct a
*b = (N - Nfix)*(*b/N); // Correct b
}

void ChangeWeight(RooRealVar *Nmix, RooRealVar *Ngas, RooRealVar *Nbk, double wmix,double c, int N){
Nmix->setVal(N*wmix); // Nmix is changing, N*wmix
Ngas->setVal(N - wmix*N - c*N); // Ngas are the remaining events, after subtracting the cosmic rate c*N
Nbk->setVal(N*c); // Fix the background 
}

void PrintInfo(RooRealVar Nmix_f,RooRealVar Nuw_f,RooRealVar Nbk_f,RooRealVar Nmix_t,RooRealVar Nuw_a,RooRealVar Nbk_a, double wmix, int i){
std::cout << "\nEvent LOOP N°" << i  <<std::endl;
	std::cout << "weight Mix: " << wmix << std::endl;
	std::cout << "fit:      " << "Nmix: " << Nmix_f.getVal() << " Nuw: " << Nuw_f.getVal() << " Nbk " << Nbk_f.getVal() << std::endl;
	std::cout << "Expected: " << "Nmix: " << Nmix_t.getVal() << " Nuw: " << Nuw_a.getVal() << " Nbk: " << Nbk_a.getVal() << std::endl;
	}

#endif
