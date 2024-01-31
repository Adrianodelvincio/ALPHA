#ifndef ANALYSISLINESHAPE_H
#define ANALYSISLINESHAPE_H
#include <iostream>
#include <string>
#include <fstream>
#include "TString.h"
#include "TSystem.h"

std::vector<std::string> getFiles(TString);
std::vector<std::string> getFiles(int start, int stop, TString Directory){
	std::vector<string> elenco;
	TString Namefile1 = "LoopDataPdf1_";
	TString Namefile2 = "LoopDataPdf2_";
	//std::cout << "Creating list of files" << std::endl;
	
	for(int i = start; i <= stop; i++){
		TString endfile = TString::Format("%d.root", i);
		TString namefile1 = Directory + Namefile1 + endfile;
		TString namefile2 = Directory + Namefile2 + endfile;
		/*
		if(i%10 == 0){
		std::cout << namefile1 << std::endl;
		std::cout << namefile2 << std::endl;
		}
		*/
		if(!gSystem->AccessPathName(namefile1) && !gSystem->AccessPathName(namefile2)){
			std::string lastfile1(namefile1.Data());
			std::string lastfile2(namefile2.Data());
			elenco.push_back(lastfile1);
			elenco.push_back(lastfile2);
		}else { break;}
	}
	return elenco;
}

bool CheckEmpty(TString RunName){
	std::string name((RunName).Data());
	ifstream stream(name); 
	string firstLine,str; 
	getline(stream, firstLine);
	if(!std::getline(stream,str)){
		std::cout << "File EMPTY" << std::endl;
		return true;
	}
	else{
		return false;
	}
}

double mean(const std::vector<double> v){
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	return sum/ v.size();
}

double stdev(const std::vector<double> v){
	double accum = 0.0;
	double m = mean(v);
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
    	accum += (d - m) * (d - m);
	});
	return sqrt(accum/(v.size() - 1));
}

double ResSquare(const std::vector<double> v){
	double result = 0;
	for(int i = 0; i < v.size(); i++){
		result += v[i]*v[i]; 
	}
	return result/v.size();
}

double Corr(const std::vector<double> v1, const std::vector<double> v2){
	double m1 = mean(v1);
	double m2 = mean(v2);
	double sigma0 = stdev(v1);
	double sigma1 = stdev(v2);
	double accum = 0.0;
	for(int i =0; i < v1.size(); i++){
		accum += (v1[i] - m1)*(v2[i] - m2);
	}
	if(sigma0 != 0.0 && sigma1 != 0.0){
	accum = accum/(sigma0*sigma1);}
	else{
	accum = 0;
	}
	return accum/v1.size();
}

#endif
