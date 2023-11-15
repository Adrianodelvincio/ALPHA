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
	
	auto cdf = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r68949_cosmics.vertex.csv");
	auto cdf1 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69177_cosmics.vertex.csv");
	auto cdf2 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69207_cosmics.vertex.csv");
	auto cdf3 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69219_cosmics.vertex.csv");
	
	cdf.Snapshot("myTree", "DataSetROOT/r68949_cosmics.vertex.root");
	cdf1.Snapshot("myTree", "DataSetROOT/r69177_cosmics.vertex.root");
	cdf2.Snapshot("myTree", "DataSetROOT/r69207_cosmics.vertex.root");
	cdf3.Snapshot("myTree", "DataSetROOT/r69219_cosmics.vertexroot");
	// REAL DATA
	std::cout << "Real Data: " << std::endl;
	TString front_file = "Dataset/";
	TString RealData[3] = {"r68465", "r68481", "r68489"};
	numFile = 1;
	for (int i = 0; i < 3 ; i ++){
		numFile = 1;
		while(numFile < 100){
		TString file_end = TString::Format("_uw_exp_freq%d.vertex.csv", numFile);
		std::cout << "Processing file: " << front_file + RealData[i] + file_end << std::endl;
			if(!gSystem->AccessPathName(front_file + RealData[i] + file_end)){
			std::cout << "File found, proceding to the conversion" << std::endl;
			}
			else{
			std::cout << "File NOT found, loop ends" << std::endl;
			break;
			}
		numFile += 1;
		};
	};
	/*
	ROOT::RDataFrame total_rdf("myTree", {"DataSetROOT/r68465_f1.root", "DataSetROOT/r68465_f2.root", "DataSetROOT/r68465_f3.root", "DataSetROOT/r68465_f4.root", "DataSetROOT/r68465_f5.root"});
	total_rdf.Snapshot("myTree", "DataSetROOT/r68465_cut1.root");

	
	ROOT::RDataFrame ttotal_rdf("myTree", {"DataSetROOT/r68481_f1.root","DataSetROOT/r68481_f2.root","DataSetROOT/r68481_f3.root","DataSetROOT/r68481_f4.root","DataSetROOT/r68481_f5.root","DataSetROOT/r68481_f6.root","DataSetROOT/r68481_f7.root","DataSetROOT/r68481_f8.root"});
	ttotal_rdf.Snapshot("myTree", "DataSetROOT/r68481_cut1.root");
	
	
	ROOT::RDataFrame tttotal_rdf("myTree", {"DataSetROOT/r68489_f1.root","DataSetROOT/r68489_f2.root","DataSetROOT/r68489_f3.root","DataSetROOT/r68489_f4.root","DataSetROOT/r68489_f5.root","DataSetROOT/r68489_f6.root","DataSetROOT/r68489_f7.root","DataSetROOT/r68489_f8.root","DataSetROOT/r68489_f9.root"});
	tttotal_rdf.Snapshot("myTree", "DataSetROOT/r68489_cut1.root");
	
	
	ROOT::RDataFrame ttttotal_rdf("myTree", {"DataSetROOT/r68498_f1.root","DataSetROOT/r68498_f2.root","DataSetROOT/r68498_f3.root","DataSetROOT/r68498_f4.root","DataSetROOT/r68498_f5.root","DataSetROOT/r68498_f6.root","DataSetROOT/r68498_f7.root"});
	ttttotal_rdf.Snapshot("myTree", "DataSetROOT/r68498_cut1.root");
	*/
		
}
