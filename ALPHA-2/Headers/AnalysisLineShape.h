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
	std::cout << "Creating list of files" << std::endl;
	
	for(int i = start; i <= stop; i++){
		TString endfile = TString::Format("%d.root", i);
		TString namefile1 = Directory + Namefile1 + endfile;
		TString namefile2 = Directory + Namefile2 + endfile;
		if(i%10 == 0){
		std::cout << namefile1 << std::endl;
		std::cout << namefile2 << std::endl;
		}
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
#endif
