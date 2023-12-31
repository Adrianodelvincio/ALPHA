#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include <ROOT/RDataFrame.hxx>


double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);
double firstWithVeto(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold);
double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf);
double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf);
double constFrac(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction);
double sumNeighbors(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold, int N);

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		if(histpdf->GetBinContent(i) > threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
		onset = histpdf->GetBinCenter(i);
	}
	std::cout << "Threshold: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double firstWithVeto(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) >= threshold){
			if(histpdf->GetBinContent(i+1) >= threshold){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	std::cout << "Veto: " << " i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}



double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) > 0){
			if(histpdf->GetBinContent(i+1) > 1){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	std::cout << "2017: " << " i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){
	// bin 0 is underflow
	double onset = 0;			// onset value 
	double bin = histpdf->GetNbinsX();	// bin onset
	//std::cout << "last bin " << bin << " content: " << histpdf->GetBinContent(bin) << std::endl;
	//std::cout << "Loop" << std::endl;
	for(int i = histpdf->GetNbinsX(); i >= 1; --i){
		//std::cout << "current bin: " << i << " content: " << histpdf->GetBinContent(i) << std::endl;
		//std::cout << "   following bin: " << i -1 << " content: " <<  histpdf->GetBinContent(i -1) << std::endl;
		if(histpdf->GetBinContent(i) < 3 && histpdf->GetBinContent(i-1) < 2){
			bin = i;
			onset = histpdf->GetBinCenter(i);
			break;
		}
		if(i == 1){
		onset = histpdf->GetBinCenter(i);
		bin = 1;
		}
	}
	std::cout << "2017 reversed" << " i: " << bin << " frequency: " << histpdf->GetBinCenter(bin)  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double constFrac(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction){
// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	double threshold = fraction*histpdf->GetMaximum();
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		if(histpdf->GetBinContent(i) > threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
	}
	std::cout << "ConstFrac: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;

}


double sumNeighbors(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold, int N){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	threshold =  3*threshold + N*sqrt(3)*sqrt(threshold);
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		double sum = histpdf->GetBinContent(i) + histpdf->GetBinContent(i+1) + histpdf->GetBinContent(i+2);
		if(sum > threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
		onset = histpdf->GetBinCenter(i);
	}
	std::cout << "sumNeighbors: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}
#endif
