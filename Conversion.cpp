#include <iostream>


void Conversion(){
	
	// Mixing
	auto list1 = {"Spectroscopy/Control/intrap/r68814_mixing.vertex.csv", "Spectroscopy/Control/intrap/r68839_mixing.vertex.csv", "Spectroscopy/Control/intrap/r68859_mixing.vertex.csv", "Spectroscopy/Control/intrap/r68871_mixing.vertex.csv", "Spectroscopy/Control/intrap/r68903_mixing.vertex.csv", "Spectroscopy/Control/intrap/r68905_mixing.vertex.csv", "Spectroscopy/Control/intrap/r68927_mixing.vertex.csv", "Spectroscopy/Control/intrap/r69126_mixing.vertex.csv","Spectroscopy/Control/intrap/r69142_mixing.vertex.csv" };
	
	auto rdf = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68814_mixing.vertex.csv");
	auto rdf1 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68839_mixing.vertex.csv");
	auto rdf2 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68859_mixing.vertex.csv");
	auto rdf3 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68871_mixing.vertex.csv");
	auto rdf4 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68903_mixing.vertex.csv");
	auto rdf5 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68905_mixing.vertex.csv");
	auto rdf6 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68927_mixing.vertex.csv");
	auto rdf7 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r69126_mixing.vertex.csv");
	auto rdf8 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r69142_mixing.vertex.csv");
	
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
	//"Spectroscopy/Control/cosmic/" "_cosmics.vertex.csv" r68949 r69177 r69207 r69219
	 
	auto cdf = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/cosmic/r68949_cosmics.vertex.csv");
	auto cdf1 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/cosmic/r69177_cosmics.vertex.csv");
	auto cdf2 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/cosmic/r69207_cosmics.vertex.csv");
	auto cdf3 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/cosmic/r69219_cosmics.vertex.csv");
	
	cdf.Snapshot("myTree", "DataSetROOT/r68949_cosmics.vertex.root");
	cdf1.Snapshot("myTree", "DataSetROOT/r69177_cosmics.vertex.root");
	cdf2.Snapshot("myTree", "DataSetROOT/r69207_cosmics.vertex.root");
	cdf3.Snapshot("myTree", "DataSetROOT/r69219_cosmics.vertexroot");
	
	// UWlosses
	
	auto udf = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68814_uwlosses_160.vertex.csv");
	auto udf1 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68839_uwlosses_160.vertex.csv");
	auto udf2 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68859_uwlosses_160.vertex.csv");
	auto udf3 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68871_uwlosses_160.vertex.csv");
	auto udf4 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68903_uwlosses_160.vertex.csv");
	auto udf5 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68905_uwlosses_160.vertex.csv");
	auto udf6 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r68927_uwlosses_160.vertex.csv");
	auto udf7 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r69126_uwlosses_160.vertex.csv");
	auto udf8 = ROOT::RDF::MakeCsvDataFrame("Spectroscopy/Control/intrap/r69142_uwlosses_160.vertex.csv");
	
	udf.Snapshot("myTree", "DataSetROOT/r68814_uwlosses_160.vertex.root");
	udf1.Snapshot("myTree", "DataSetROOT/r68839_uwlosses_160.vertex.root");
	udf2.Snapshot("myTree", "DataSetROOT/r68859_uwlosses_160.vertex.root");
	udf3.Snapshot("myTree", "DataSetROOT/r68871_uwlosses_160.vertex.root");
	udf4.Snapshot("myTree", "DataSetROOT/r68903_uwlosses_160.vertex.root");
	udf5.Snapshot("myTree", "DataSetROOT/r68905_uwlosses_160.vertex.root");
	udf6.Snapshot("myTree", "DataSetROOT/r68927_uwlosses_160.vertex.root");
	udf7.Snapshot("myTree", "DataSetROOT/r69126_uwlosses_160.vertex.root");
	udf8.Snapshot("myTree", "DataSetROOT/r69142_uwlosses_160.vertex.root");
	
}
