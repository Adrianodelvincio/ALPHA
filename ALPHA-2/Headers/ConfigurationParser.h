#ifndef CONFIGURATIONPARSER_H
#define CONFIGURATIONPARSER_H
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
					// std::cout << line << std::endl; // Print Parsed Line
					if(line.find("#")!=std::string::npos) continue;
					char sBuffer[20]; float buffer;
					if (sscanf(line.c_str(), "%s = %f", sBuffer, &buffer) == 2) {
						if (strcmp(sBuffer, "pwall_cb") == 0) 		{ pwall_cb = buffer;}
						if (strcmp(sBuffer, "pwall_ad") == 0) 		{ pwall_ad = buffer;}
						if (strcmp(sBuffer, "C") == 0) 				{ C = buffer;}
						if (strcmp(sBuffer, "x_cb_start") == 0) 	{ x_cb_start = buffer;}
						if (strcmp(sBuffer, "x_cb_end") == 0) 		{ x_cb_end = buffer;}
						if (strcmp(sBuffer, "x_cb_peak") == 0) 		{ x_cb_peak = buffer;}
						if (strcmp(sBuffer, "x_da_start") == 0) 	{ x_da_start = buffer;}
						if (strcmp(sBuffer, "x_da_peak") == 0) 		{ x_da_peak = buffer;}
						if (strcmp(sBuffer, "x_da_end") == 0) 		{ x_da_end = buffer;}
						if (strcmp(sBuffer, "CosmicRate") == 0) 	{ CosmicRate = buffer;}
						if (strcmp(sBuffer, "Efficiency") == 0) 	{ Efficiency = buffer;}
						if (strcmp(sBuffer, "delay") == 0) 			{ delay = buffer;}
						if (strcmp(sBuffer, "FrequencyStep") == 0) 	{ FrequencyStep = buffer;}
					} // fill the parameters
					int IntBuffer;
					if (sscanf(line.c_str(), "%s = %d", sBuffer, &IntBuffer) == 2) {
						if (strcmp(sBuffer, "Nstack") == 0) 		{ Nstack = IntBuffer;}
						if (strcmp(sBuffer, "NHbar") == 0) 			{ NHbar = IntBuffer;}
						if (strcmp(sBuffer, "Repetition") == 0) 	{ Repetition = IntBuffer;}
						if (strcmp(sBuffer, "TimeStep") == 0) 		{ TimeStep = IntBuffer;}
						if (strcmp(sBuffer, "SweepStep") == 0) 		{ SweepStep = IntBuffer;}
					}
				} // End read line loop
		} // End Constructor
		
		void Print(){
			std::cout << "\n" << std::endl;
			std::cout << "Nstack        " << Nstack << std::endl;
			std::cout << "NHbar         " << NHbar<< std::endl;
			std::cout << "Repetition    " << Repetition  << std::endl;
			std::cout << "TimeStep      " << TimeStep  << std::endl;
			std::cout << "SweepStep     " << SweepStep<< std::endl;
			std::cout << "FrequencyStep " << FrequencyStep << std::endl;
			std::cout << "Efficiency    " << Efficiency << std::endl;
			std::cout << "CosmicRate    " << CosmicRate  << std::endl;
			std::cout << "x_cb_start    " << x_cb_start << std::endl;
			std::cout << "x_cb_end      " << x_cb_end << std::endl;
			std::cout << "x_cb_peak     " << x_cb_peak << std::endl;
			std::cout << "x_da_start    " << x_da_start << std::endl;
			std::cout << "x_da_peak     " << x_da_peak << std::endl;
			std::cout << "x_da_end      " << x_da_end << std::endl;
			std::cout << "C             " << C << std::endl;
			std::cout << "pwall_ad      " << pwall_ad << std::endl;
			std::cout << "pwall_cb      " <<  pwall_cb << std::endl;
			std::cout << "delay         " <<  delay << std::endl;
			std::cout << "\n" << std::endl;
		}
	int Nstack;
	int NHbar;
	int Repetition;
	int TimeStep;
	int SweepStep;
	double FrequencyStep;
	double x_cb_start;
	double x_cb_end;
	double x_cb_peak;
	double x_da_start;
	double x_da_peak;
	double x_da_end;
	double CosmicRate;
	double Efficiency;
	double pwall_cb;
	double pwall_ad;
	double C;
	double delay;
};
#endif
