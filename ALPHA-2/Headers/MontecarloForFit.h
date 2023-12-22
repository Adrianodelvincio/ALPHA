#ifndef MONTECARLOFORFIT_H
#define MONTECARLOFORFIT_H
#include <iostream>
#include "RooRealVar.h"
#include <numeric>
#include <fstream>
#include <iostream>
#include "TString.h"


class ReadConfFile{
	public:
		ReadConfFile();
		ReadConfFile(TString fileName){
			std::string file(fileName.Data());
			std::ifstream myfile (file.c_str());	// Open the file
			std::string line; 						// Parse each line
				
			while(myfile.is_open() && getline(myfile,line)){
				// If the line is a comment
				if(line.find("#")!=std::string::npos) continue;
				char sBuffer[20]; float buffer;
				if (sscanf(line.c_str(), "%s = %f", sBuffer, &buffer) == 2) {
					if (strcmp(sBuffer, "a") == 0) 	{ a = buffer;}
					if (strcmp(sBuffer, "b") == 0) 	{ b = buffer;}
					if (strcmp(sBuffer, "Ncosmic") == 0) 	{ Ncosmic = buffer;}
					if (strcmp(sBuffer, "mu") == 0) 	{ mu = buffer;}
					if (strcmp(sBuffer, "sigWall") == 0) 	{ sigWall = buffer;}
					if (strcmp(sBuffer, "sigRay") == 0) 	{ sigRay = buffer;}
				} // fill the parameters
				int IntBuffer;
				if (sscanf(line.c_str(), "%s = %d", sBuffer, &IntBuffer) == 2) {
					if (strcmp(sBuffer, "Nloop") == 0) 	{ Nloop = IntBuffer;}
					if (strcmp(sBuffer, "N") == 0) 		{ N = IntBuffer;}
				}
			} // End read line loop
		} // End Constructor
		
		void Print(){
			std::cout << "\n" << std::endl;
			std::cout << "a:     " << a	<< std::endl;
			std::cout << "b:     " << b	<< std::endl;
			std::cout << "N:     " << N	<< std::endl;
			std::cout << "Nloop: " << Nloop	<< std::endl;
			std::cout << "Ncosmic: " << Ncosmic << std::endl;
			std::cout << "\n" << std::endl;
		}
		int Nloop,N;
		double a,b,Ncosmic,sigWall,mu,sigRay;
};



void SetVectors(int i,std::vector<double> &mix,std::vector<double> &uw, std::vector<double> &bk, RooRealVar a,RooRealVar b,RooRealVar c,RooRealVar aa,RooRealVar bb,RooRealVar cc){
		mix.push_back(a.getVal() - aa.getVal());
		uw.push_back(b.getVal() - bb.getVal());
		bk.push_back(c.getVal() - cc.getVal());}
		

Double_t media(std::vector<Double_t>&v) {
if(v.size()<=0.) return 0.;
return (Double_t) std::accumulate(v.begin(), v.end(), 0.0) / (Double_t)v.size();
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
	std::cout << "Expected: " << "Nmix: " << Nmix_t.getVal() << " Nuw: " << Nuw_a.getVal() << " Nbk: " << Nbk_a.getVal() << "\n" <<std::endl;
	}

#endif
