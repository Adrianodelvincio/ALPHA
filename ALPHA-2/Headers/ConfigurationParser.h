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
						if (strcmp(sBuffer, "WallComponent_cb") == 0) 		{ WallComponent_cb = buffer;}
						if (strcmp(sBuffer, "WallComponent_ad") == 0) 		{ WallComponent_ad = buffer;}
						if (strcmp(sBuffer, "C") == 0) 			{ C = buffer;}
						if (strcmp(sBuffer, "x_cb_start") == 0) 	{ x_cb_start = buffer;}
						if (strcmp(sBuffer, "x_cb_end") == 0) 		{ x_cb_end = buffer;}
						if (strcmp(sBuffer, "x_cb_peak") == 0) 		{ x_cb_peak = buffer;}
						if (strcmp(sBuffer, "x_da_start") == 0) 	{ x_da_start = buffer;}
						if (strcmp(sBuffer, "x_da_peak") == 0) 		{ x_da_peak = buffer;}
						if (strcmp(sBuffer, "x_da_end") == 0) 		{ x_da_end = buffer;}
						if (strcmp(sBuffer, "CosmicRate") == 0) 	{ CosmicRate = buffer;}
						if (strcmp(sBuffer, "Efficiency") == 0) 	{ Efficiency = buffer;}
						if (strcmp(sBuffer, "delay") == 0) 		{ delay = buffer;}
						if (strcmp(sBuffer, "FrequencyStep") == 0) 	{ FrequencyStep = buffer;}
						if (strcmp(sBuffer, "Norm_cb") == 0)		{ Norm_cb = buffer;}
						if (strcmp(sBuffer, "x0_cb") == 0)		{ x0_cb = buffer;}
						if (strcmp(sBuffer, "sigma0_cb") == 0)		{ sigma0_cb = buffer;}
						if (strcmp(sBuffer, "sigma1_cb") == 0)		{ sigma1_cb = buffer;}
						if (strcmp(sBuffer, "k0_cb") == 0)		{ k0_cb = buffer;}
						if (strcmp(sBuffer, "k1_cb") == 0)		{ k1_cb = buffer;}
						if (strcmp(sBuffer, "x0_da") == 0)		{ x0_da = buffer;}
						if (strcmp(sBuffer, "sigma0_da") == 0)		{ sigma0_da = buffer;}
						if (strcmp(sBuffer, "sigma1_da") == 0)		{ sigma1_da = buffer;}
						if (strcmp(sBuffer, "k0_da") == 0)		{ k0_da = buffer;}
						if (strcmp(sBuffer, "k1_da") == 0)		{ k1_da = buffer;}
						if (strcmp(sBuffer, "NHbar") == 0) 		{ NHbar = buffer;}
						if (strcmp(sBuffer, "Norm_da") == 0)		{ Norm_da = buffer;}
						if (strcmp(sBuffer, "delta_right") == 0)	{ delta_right = buffer;}
						if (strcmp(sBuffer, "delta_left") == 0)		{ delta_left = buffer;}
						if (strcmp(sBuffer, "FrequencyMaxBeforeOnset") == 0){ FrequencyMaxBeforeOnset = buffer;}
						if (strcmp(sBuffer, "FrequencyMinBeforeOnset") == 0){ FrequencyMinBeforeOnset = buffer;}
						if (strcmp(sBuffer, "Bdrift") == 0)		{ Bdrift = buffer;}
						if (strcmp(sBuffer, "sigmaDrift") == 0)		{ sigmaDrift = buffer;}
						
					} // fill the parameters
					int IntBuffer;
					if (sscanf(line.c_str(), "%s = %d", sBuffer, &IntBuffer) == 2) {
						if (strcmp(sBuffer, "Nstack") == 0) 		{ Nstack = IntBuffer;}
						if (strcmp(sBuffer, "Repetition") == 0) 	{ Repetition = IntBuffer;}
						if (strcmp(sBuffer, "TimeStep") == 0) 		{ TimeStep = IntBuffer;}
						if (strcmp(sBuffer, "SweepStep") == 0) 		{ SweepStep = IntBuffer;}
						if (strcmp(sBuffer, "BinBeforeOnset") == 0) 	{ BinBeforeOnset = IntBuffer;}
						if (strcmp(sBuffer, "TotalStep") == 0) 		{ TotalStep = IntBuffer;}
					}
				} // End read line loop
		} // End Constructor
		
		void Print(){
			std::cout << "\n" << std::endl;
			std::cout << "Nstack        " << Nstack << std::endl;
			std::cout << "NHbar         " << NHbar<< std::endl;
			std::cout << "Repetition    " << Repetition  << std::endl;
			std::cout << "TimeStep      " << TimeStep  << std::endl;
			std::cout << "SweepStep     " << SweepStep << std::endl;
			std::cout << "TotalStep     " << TotalStep << std::endl;
			std::cout << "BinBefore     " << BinBeforeOnset << std::endl;
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
			std::cout << "pwall_ad      " << WallComponent_ad << std::endl;
			std::cout << "pwall_cb      " << WallComponent_cb << std::endl;
			std::cout << "delay         " << delay << std::endl;
			std::cout << "Bdrift        " << Bdrift << std::endl;
			
			std::cout << "Cruijff Paramters		" << std::endl;
			std::cout << "x0_cb     " << x0_cb << std::endl;
			std::cout << "sigma0_cb	" << sigma0_cb << std::endl;
			std::cout << "sigma1_cbì " << sigma1_cb << std::endl;
			std::cout << "k0_cb      " << k0_cb << std::endl;
			std::cout << "k1_cb      " << k1_cb << std::endl;
			std::cout << "Norm_cb    " << Norm_cb << std::endl;
			std::cout << "sigma0_da  " << sigma0_da << std::endl;
			std::cout << "sigma1_da  " << sigma1_da << std::endl;
			std::cout << "k0_da      " << k0_da << std::endl;
			std::cout << "k1_da      " << k1_da << std::endl;
			std::cout << "Norm_da    " << Norm_da << std::endl;
			std::cout << "\n" << std::endl;
		}
		
	// LIST OF PARAMETERS OF THE SIMULATIONS	
	int Nstack;
	int Repetition;
	int TimeStep;
	int SweepStep;
	int BinBeforeOnset;
	int TotalStep;
	double NHbar;
	double FrequencyStep;
	double x_cb_start;
	double x_cb_end;
	double x_cb_peak;
	double x_da_start;
	double x_da_peak;
	double x_da_end;
	double CosmicRate;
	double Efficiency;
	double WallComponent_cb;
	double WallComponent_ad;
	double C;
	double delay;
	double delta_left;
	double delta_right;
	double FrequencyMaxBeforeOnset;
	double FrequencyMinBeforeOnset;
	double Bdrift, sigmaDrift;
	// Parameters of the Cruijff Function
	double x0_cb;
	double sigma0_cb;
	double sigma1_cb;
	double k0_cb;
	double k1_cb;
	double Norm_cb;
	double x0_da;
	double sigma0_da;
	double sigma1_da;
	double k0_da;
	double k1_da;
	double Norm_da;
		
};
#endif
