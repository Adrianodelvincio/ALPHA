#include <iostream>


void Conversion(){
	
	// Mixing
	auto list1 = {"Control/intrap/r68814_mixing.vertex.csv", "Control/intrap/r68839_mixing.vertex.csv", "Control/intrap/r68859_mixing.vertex.csv", "Control/intrap/r68871_mixing.vertex.csv", "Control/intrap/r68903_mixing.vertex.csv", "Control/intrap/r68905_mixing.vertex.csv", "Control/intrap/r68927_mixing.vertex.csv", "Control/intrap/r69126_mixing.vertex.csv","Control/intrap/r69142_mixing.vertex.csv" };
	
	auto rdf = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68814_mixing.vertex.csv");
	auto rdf1 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68839_mixing.vertex.csv");
	auto rdf2 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68859_mixing.vertex.csv");
	auto rdf3 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68871_mixing.vertex.csv");
	auto rdf4 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68903_mixing.vertex.csv");
	auto rdf5 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68905_mixing.vertex.csv");
	auto rdf6 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68927_mixing.vertex.csv");
	auto rdf7 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r69126_mixing.vertex.csv");
	auto rdf8 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r69142_mixing.vertex.csv");
	
	rdf.Snapshot("myTree", "DataSetROOT/r68814_mixing.vertex.root");
	rdf1.Snapshot("myTree", "DataSetROOT/r68839_mixing.vertex.root");
	rdf2.Snapshot("myTree", "DataSetROOT/r68859_mixing.vertex.root");
	rdf3.Snapshot("myTree", "DataSetROOT/r68871_mixing.vertex.root");
	rdf4.Snapshot("myTree", "DataSetROOT/r68903_mixing.vertex.root");
	rdf5.Snapshot("myTree", "DataSetROOT/r68905_mixing.vertex.root");
	rdf6.Snapshot("myTree", "DataSetROOT/r68927_mixing.vertex.root");
	rdf7.Snapshot("myTree", "DataSetROOT/r69126_mixing.vertex.root");
	rdf8.Snapshot("myTree", "DataSetROOT/r69142_mixing.vertex.root");
	
	//Cosmic
	//"Control/cosmic/" "_cosmics.vertex.csv" r68949 r69177 r69207 r69219
	 
	auto cdf = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r68949_cosmics.vertex.csv");
	auto cdf1 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69177_cosmics.vertex.csv");
	auto cdf2 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69207_cosmics.vertex.csv");
	auto cdf3 = ROOT::RDF::MakeCsvDataFrame("Control/cosmic/r69219_cosmics.vertex.csv");
	
	cdf.Snapshot("myTree", "DataSetROOT/r68949_cosmics.vertex.root");
	cdf1.Snapshot("myTree", "DataSetROOT/r69177_cosmics.vertex.root");
	cdf2.Snapshot("myTree", "DataSetROOT/r69207_cosmics.vertex.root");
	cdf3.Snapshot("myTree", "DataSetROOT/r69219_cosmics.vertexroot");
	
	// UWlosses
	
	auto udf = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68814_uwlosses_160.vertex.csv");
	auto udf1 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68839_uwlosses_160.vertex.csv");
	auto udf2 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68859_uwlosses_160.vertex.csv");
	auto udf3 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68871_uwlosses_160.vertex.csv");
	auto udf4 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68903_uwlosses_160.vertex.csv");
	auto udf5 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68905_uwlosses_160.vertex.csv");
	auto udf6 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r68927_uwlosses_160.vertex.csv");
	auto udf7 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r69126_uwlosses_160.vertex.csv");
	auto udf8 = ROOT::RDF::MakeCsvDataFrame("Control/intrap/r69142_uwlosses_160.vertex.csv");
	
	udf.Snapshot("myTree", "DataSetROOT/r68814_uwlosses_160.vertex.root");
	udf1.Snapshot("myTree", "DataSetROOT/r68839_uwlosses_160.vertex.root");
	udf2.Snapshot("myTree", "DataSetROOT/r68859_uwlosses_160.vertex.root");
	udf3.Snapshot("myTree", "DataSetROOT/r68871_uwlosses_160.vertex.root");
	udf4.Snapshot("myTree", "DataSetROOT/r68903_uwlosses_160.vertex.root");
	udf5.Snapshot("myTree", "DataSetROOT/r68905_uwlosses_160.vertex.root");
	udf6.Snapshot("myTree", "DataSetROOT/r68927_uwlosses_160.vertex.root");
	udf7.Snapshot("myTree", "DataSetROOT/r69126_uwlosses_160.vertex.root");
	udf8.Snapshot("myTree", "DataSetROOT/r69142_uwlosses_160.vertex.root");
	
	
	// REAL DATA
	auto realf = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq1.vertex.csv");
	auto realf1 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq1.vertex.csv");
	auto realf2 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq2.vertex.csv");
	auto realf3 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq3.vertex.csv");
	auto realf4 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq4.vertex.csv");
	auto realf5 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq5.vertex.csv");
	realf1.Snapshot("myTree", "DataSetROOT/r68465_f1.root");
	realf2.Snapshot("myTree", "DataSetROOT/r68465_f2.root");
	realf3.Snapshot("myTree", "DataSetROOT/r68465_f3.root");
	realf4.Snapshot("myTree", "DataSetROOT/r68465_f4.root");
	realf5.Snapshot("myTree", "DataSetROOT/r68465_f5.root");
	
	ROOT::RDataFrame total_rdf("myTree", {"DataSetROOT/r68465_f1.root", "DataSetROOT/r68465_f2.root", "DataSetROOT/r68465_f3.root", "DataSetROOT/r68465_f4.root", "DataSetROOT/r68465_f5.root"});
	total_rdf.Snapshot("myTree", "DataSetROOT/r68465_cut1.root");

	auto realff1 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq1.vertex.csv");
	auto realff2 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq2.vertex.csv");
	auto realff3 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq3.vertex.csv");
	auto realff4 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq4.vertex.csv");
	auto realff5 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq5.vertex.csv");
	auto realff6 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq6.vertex.csv");
	auto realff7 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq7.vertex.csv");
	auto realff8 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68481_uw_exp_freq8.vertex.csv");
	realff1.Snapshot("myTree", "DataSetROOT/r68481_f1.root");
	realff2.Snapshot("myTree", "DataSetROOT/r68481_f2.root");
	realff3.Snapshot("myTree", "DataSetROOT/r68481_f3.root");
	realff4.Snapshot("myTree", "DataSetROOT/r68481_f4.root");
	realff5.Snapshot("myTree", "DataSetROOT/r68481_f5.root");
	realff6.Snapshot("myTree", "DataSetROOT/r68481_f6.root");
	realff7.Snapshot("myTree", "DataSetROOT/r68481_f7.root");
	realff8.Snapshot("myTree", "DataSetROOT/r68481_f8.root");
	
	ROOT::RDataFrame ttotal_rdf("myTree", {"DataSetROOT/r68481_f1.root","DataSetROOT/r68481_f2.root","DataSetROOT/r68481_f3.root","DataSetROOT/r68481_f4.root","DataSetROOT/r68481_f5.root","DataSetROOT/r68481_f6.root","DataSetROOT/r68481_f7.root","DataSetROOT/r68481_f8.root"});
	ttotal_rdf.Snapshot("myTree", "DataSetROOT/r68481_cut1.root");
	
	
	auto realfff1 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq1.vertex.csv");
	auto realfff2 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq2.vertex.csv");
	auto realfff3 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq3.vertex.csv");
	auto realfff4 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq4.vertex.csv");
	auto realfff5 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq5.vertex.csv");
	auto realfff6 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq6.vertex.csv");
	auto realfff7 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq7.vertex.csv");
	auto realfff8 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq8.vertex.csv");
	auto realfff9 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68489_uw_exp_freq9.vertex.csv");
	realfff1.Snapshot("myTree", "DataSetROOT/r68489_f1.root");
	realfff2.Snapshot("myTree", "DataSetROOT/r68489_f2.root");
	realfff3.Snapshot("myTree", "DataSetROOT/r68489_f3.root");
	realfff4.Snapshot("myTree", "DataSetROOT/r68489_f4.root");
	realfff5.Snapshot("myTree", "DataSetROOT/r68489_f5.root");
	realfff6.Snapshot("myTree", "DataSetROOT/r68489_f6.root");
	realfff7.Snapshot("myTree", "DataSetROOT/r68489_f7.root");
	realfff8.Snapshot("myTree", "DataSetROOT/r68489_f8.root");
	realfff8.Snapshot("myTree", "DataSetROOT/r68489_f9.root");
	
	ROOT::RDataFrame tttotal_rdf("myTree", {"DataSetROOT/r68489_f1.root","DataSetROOT/r68489_f2.root","DataSetROOT/r68489_f3.root","DataSetROOT/r68489_f4.root","DataSetROOT/r68489_f5.root","DataSetROOT/r68489_f6.root","DataSetROOT/r68489_f7.root","DataSetROOT/r68489_f8.root","DataSetROOT/r68489_f9.root"});
	tttotal_rdf.Snapshot("myTree", "DataSetROOT/r68489_cut1.root");
	
	auto realffff1 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq1.vertex.csv");
	auto realffff2 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq2.vertex.csv");
	auto realffff3 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq3.vertex.csv");
	auto realffff4 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq4.vertex.csv");
	auto realffff5 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq5.vertex.csv");
	auto realffff6 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq6.vertex.csv");
	auto realffff7 = ROOT::RDF::MakeCsvDataFrame("Dataset/r68498_uw_exp_freq7.vertex.csv");
	
	realffff1.Snapshot("myTree", "DataSetROOT/r68498_f1.root");
	realffff2.Snapshot("myTree", "DataSetROOT/r68498_f2.root");
	realffff3.Snapshot("myTree", "DataSetROOT/r68498_f3.root");
	realffff4.Snapshot("myTree", "DataSetROOT/r68498_f4.root");
	realffff5.Snapshot("myTree", "DataSetROOT/r68498_f5.root");
	realffff6.Snapshot("myTree", "DataSetROOT/r68498_f6.root");
	realffff7.Snapshot("myTree", "DataSetROOT/r68498_f7.root");
	
	
	ROOT::RDataFrame ttttotal_rdf("myTree", {"DataSetROOT/r68498_f1.root","DataSetROOT/r68498_f2.root","DataSetROOT/r68498_f3.root","DataSetROOT/r68498_f4.root","DataSetROOT/r68498_f5.root","DataSetROOT/r68498_f6.root","DataSetROOT/r68498_f7.root"});
	ttttotal_rdf.Snapshot("myTree", "DataSetROOT/r68498_cut1.root");
	
		
}
