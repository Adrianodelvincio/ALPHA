#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include <ROOT/RDataFrame.hxx>

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold, double background){ // WITH COSMIC BACKGROUND
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		if(histpdf->GetBinContent(i) > round(threshold + background)){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
		onset = histpdf->GetBinCenter(i);
	}
	//std::cout << "Threshold: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){ // WITHOUT COSMIC BACKGROUND
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
	//std::cout << "Threshold: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
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
	//std::cout << "Veto: " << " i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background){ //WITH COSMIC BACKGROUND
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) > round(0 + background)){
			if(histpdf->GetBinContent(i+1) > round(1 + background)){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	//std::cout << "2017: " << " i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){ //WITHOUT COSMIC BACKGROUND
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 0; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) > 0){
			if(histpdf->GetBinContent(i+1) > 1){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
	}
	std::cout << "2017: " << " frequence " << histpdf->GetBinCenter(bin)  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background){ // WITH COSMIC BACKGROUND
	// bin 0 is underflow
	double onset = 0;			// onset value 
	double bin = histpdf->GetNbinsX();	// bin onset
	for(int i = histpdf->GetMaximumBin(); i >= 1; --i){
		if(histpdf->GetBinContent(i) < round(3 + background) && histpdf->GetBinContent(i-1) < round(2 + background)){
			bin = i;
			onset = histpdf->GetBinCenter(i);
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	//std::cout << "2017 reversed" << " i: " << bin << " frequency: " << histpdf->GetBinCenter(bin)  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){ // WITHOUT COSMIC BACKGROUND
	// bin 0 is underflow
	double onset = 0;			// onset value 
	double bin = histpdf->GetNbinsX();	// bin onset
	for(int i = histpdf->GetMaximumBin(); i >= 1; --i){
		if(histpdf->GetBinContent(i) < 3 && histpdf->GetBinContent(i-1) < 2){
			bin = i;
			onset = histpdf->GetBinCenter(i);
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	//std::cout << "2017 reversed" << " i: " << bin << " frequency: " << histpdf->GetBinCenter(bin)  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double constFrac(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction, double background){ // WITH BACKGROUND
// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	double threshold = fraction*(histpdf->GetMaximum() - background);
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		if(histpdf->GetBinContent(i) > (threshold + background)){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
	}
	//std::cout << "ConstFrac: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;

}

double constFrac(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction){ // WITHOUT BACKGROUND
// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	double threshold = fraction*histpdf->GetMaximum();
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		bin = i;
		if((histpdf->GetBinContent(i)) > threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
	}
	//std::cout << "ConstFrac: " << threshold << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double sumNeighbors(ROOT::RDF::RResultPtr<TH1D> histpdf, int N, double background){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 2;		// bin onset
	double threshold;
	threshold =  3*background + N*sqrt(3)*sqrt(background);
	for(int i = 2; i < histpdf->GetNbinsX() - 1; i++){
		bin = i;
		double sum = histpdf->GetBinContent(i -1) + histpdf->GetBinContent(i) + histpdf->GetBinContent(i+1);
		if(sum > threshold){
			onset = histpdf->GetBinCenter(i);
			bin = i;
			break;
		}
		onset = histpdf->GetBinCenter(i);
	}
	//std::cout << "sumNeighbors: " << threshold << " background: " << background << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double runningDiff(ROOT::RDF::RResultPtr<TH1D> histpdf){
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	double threshold = 2;
	for(int i = 1; i < histpdf->GetNbinsX() - 2; i++){
		double diff1 = histpdf->GetBinContent(i + 1) - histpdf->GetBinContent(i);
		double diff2 = histpdf->GetBinContent(i + 2) - histpdf->GetBinContent(i+1);
		if(diff1 >= threshold && diff2 >= threshold){
			onset = histpdf->GetBinCenter(i+1);
			bin = i +1;
			break;
		}
		onset = histpdf->GetBinCenter(i);
	}
	std::cout << "runningDiff: " <<  " frequency " << histpdf->GetBinCenter(bin) << std::endl;
	return onset;
	

}

double hybrid_cfSum(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction ,double background){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	double threshold;
	threshold =  3*(background + ( histpdf->GetMaximum() - background )*fraction);
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
	//std::cout << "hybrid: " << threshold << " background: " << background << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;	
}
#endif
