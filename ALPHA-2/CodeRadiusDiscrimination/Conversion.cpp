#include <iostream>


void Conversion(){
	// Generazione file .root Template
	TString Mixing_Files[9] = {"r68814", "r68839", "r68859", "r68871", "r68903", "r68905", "r68927", "r69126","r69142" };
	TString folderCsv = "Control/intrap/";
	TString endnameMix = "_mixing.vertex.csv";
	TString endnameUw = "_uwlosses_160.vertex.csv";
	TString folder = TString::Format("DataSetROOT/");
	TString endfilerootMix = "_mixing.vertex.root";
	TString endfilerootUw = "_uwlosses_160.vertex.root";
	
	int numFile = 0;
	while(numFile <= 8){
		// Mixing
		std::cout << "Processing file : " << folderCsv+Mixing_Files[numFile]+endnameMix << std::endl;
		if(!gSystem->AccessPathName(folderCsv+Mixing_Files[numFile]+endnameMix)){
			std::cout << "File found, creating .root file" << std::endl;
			if(gSystem->AccessPathName(folder + Mixing_Files[numFile] + endfilerootMix)){
			auto rdf = ROOT::RDF::MakeCsvDataFrame(folderCsv + Mixing_Files[numFile] + endnameMix);
			rdf.Snapshot("myTree", folder + Mixing_Files[numFile] + endfilerootMix);}
			else{
			std::cout << "!!!File .root already existing!!!" << std::endl;
			}
		}
		// Residual Gas
		std::cout << "Processing file : " << folderCsv+Mixing_Files[numFile]+endnameUw << std::endl;
		if(!gSystem->AccessPathName(folderCsv+Mixing_Files[numFile]+endnameUw)){
			std::cout << "File found, creating .root file" << std::endl;
			if(gSystem->AccessPathName(folder + Mixing_Files[numFile] + endfilerootUw)){
			auto rdf = ROOT::RDF::MakeCsvDataFrame(folderCsv + Mixing_Files[numFile] + endnameUw);
			rdf.Snapshot("myTree", folder + Mixing_Files[numFile] + endfilerootUw);}
			else{
			std::cout << "!!!File .root already existing!!!" << std::endl;
			}
		}
		
	numFile += 1;
	}
	
	// COSMICI
	auto cdf = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r68949_cosmics.vertex.csv");
	auto cdf1 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69177_cosmics.vertex.csv");
	auto cdf2 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69207_cosmics.vertex.csv");
	auto cdf3 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69219_cosmics.vertex.csv");
	if(gSystem->AccessPathName("DataSetROOT/r68949_cosmics.vertex.root")){
		cdf.Snapshot("myTree", "DataSetROOT/r68949_cosmics.vertex.root");
		cdf1.Snapshot("myTree", "DataSetROOT/r69177_cosmics.vertex.root");
		cdf2.Snapshot("myTree", "DataSetROOT/r69207_cosmics.vertex.root");
		cdf3.Snapshot("myTree", "DataSetROOT/r69219_cosmics.vertexroot");
	}
	
	// REAL DATA
	std::cout << "Real Data: " << std::endl;
	TString front_file = "Dataset/";
	TString RealData[4] = {"r68465", "r68481", "r68489", "r68498"};
	
	numFile = 1;
	for (int i = 0; i < 4 ; i ++){
		numFile = 1;
		while(numFile < 100){
		
		TString file_end = TString::Format("_uw_exp_freq%d.vertex.csv", numFile);
		std::cout << "Processing file: " << front_file + RealData[i] + file_end << std::endl;
			//Check if the file exist
			if(!gSystem->AccessPathName(front_file + RealData[i] + file_end)){	
			std::cout << "File found, proceding to the conversion" << std::endl;
			//Check Empty file
			std::string namecsv((front_file + RealData[i] + file_end).Data());
			ifstream stream(namecsv); 
			string firstLine,str; 
			getline(stream, firstLine);
			if(!std::getline(stream,str)){
			numFile += 1;
			std::cout << "File EMPTY" << std::endl;
			continue;
			} // empty files, ignore it
			// Open the Data Frame
			TString NamefileRoot = TString::Format("_f%d.root", numFile);
			auto rdf = ROOT::RDF::MakeCsvDataFrame(front_file + RealData[i] + file_end);
			rdf.Snapshot("myTree",folder + RealData[i] + NamefileRoot);
			}
			else{
			std::cout << "File NOT found, loop ends" << std::endl;
			break;
			}
		numFile += 1;
		};
	};		
}
