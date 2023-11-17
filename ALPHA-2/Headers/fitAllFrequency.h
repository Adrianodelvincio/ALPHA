#ifndef FITALLFREQUENCY_H
#define FITALLFREQUENCY_H
#include <iostream>
#include <string>
#include <fstream>
#include "TString.h"
#include "TSystem.h"
std::vector<std::string> getFiles(TString);
std::vector<std::string> getFiles(TString RunNumber){
	std::vector<string> elenco;
	TString Directory = "DataSetROOT/";
	int numFile = 1;
	std::cout << "Creating list of files" << std::endl;
	while(numFile < 100){
		TString endfile = TString::Format("_f%d.root", numFile);
		TString namefile = Directory + RunNumber + endfile;
		std::cout << namefile << std::endl;
		if(!gSystem->AccessPathName(namefile)){
			std::string lastfile(namefile.Data());
			elenco.push_back(lastfile);
		}
		else{
			break;
		}
		numFile += 1;
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
#endif
