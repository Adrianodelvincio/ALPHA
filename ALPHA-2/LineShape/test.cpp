auto FillTry(ROOT::RDataFrame &d1){
	int j = 0; 	// Variable for loop
	return d1.Define("id", [&j](){return j++;});	// Id of the events
}

void test() {

    ROOT::RDataFrame d(100); // PDF2
	auto filledFrame = FillTry(d);
	filledFrame.Snapshot("myTree", "a.root");

    // Print the content of the file, for the sake of the example
    TFile f("a.root");
    f.ls();
}
