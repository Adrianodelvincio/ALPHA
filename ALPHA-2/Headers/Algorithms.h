#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include <ROOT/RDataFrame.hxx>
#include "TH1.h"


/////////////////////////////////////////////
// OVER THRESHOLD

void FilterHistogram(TH1D* original, TH1D *filtered, int Nfilter){
	for(int i = 1; i < original->GetNbinsX() - 1; i++){
		double sum = 0;
		for(int j = 0; j < Nfilter; j++){
			if(i + j < original->GetNbinsX() - 1){
				sum += original->GetBinContent(i + j);
			}else{
				sum += 0;
			}		
		}
		filtered->SetBinContent(i, sum);
	}
}

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold, double background){ // WITH BASELINE

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
	return onset;
}

double firstOverThreshold(ROOT::RDF::RResultPtr<TH1D> histpdf, double threshold){ // WITHOUT COSMIC BACKGROUND
	
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
	return onset;
}

/////////////////////////////////////////////
// ALGORITHM 2017

// WITH BASELINE
double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background){ //WITH BASELINE
	
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
	return onset;
}

// ORIGINAL
double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){ //WITHOUT COSMIC BACKGROUND
	
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
	return onset;
}

// WITH BAKGROUND AND TWO FREE PARAMETERS
double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background, double thr1, double thr2){ //WITH BASELINE
	
	double onset = 0;	// onset value 
	double bin = 1;		// bin onset
	for(int i = 1; i < histpdf->GetNbinsX(); i++){
		if(histpdf->GetBinContent(i) > round(thr1 + background)){
			if(histpdf->GetBinContent(i+1) > round(thr2 + background)){
				onset = histpdf->GetBinCenter(i);
				bin = i;
				break;
			}
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	return onset;
}

/////////////////////////////////////////////
// REVERSED 2017

// WITH BASELINE
double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background){ // WITH BASELINE
	
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
	return onset;
}

// ORIGINAL REVERSED
double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf){ // WITHOUT COSMIC BACKGROUND
	
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
	return onset;
}

// WITH BASELINE
double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background, double Thr1, double Thr2){ // WITH BASELINE
	
	double onset = 0;			// onset value 
	double bin = histpdf->GetNbinsX();	// bin onset
	for(int i = histpdf->GetMaximumBin(); i >= 1; --i){
		if(histpdf->GetBinContent(i) < round(Thr1 + background) && histpdf->GetBinContent(i-1) < round(Thr2 + background)){
			bin = i;
			onset = histpdf->GetBinCenter(i);
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	return onset;
}

/////////////////////////////////////////////
// CONSTANT FRACTION ALGORITHM

double constFrac(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction, double background, int Nfilter){ // WITH BACKGROUND
	
	double onset = 0;   // value of the frequency onset
	double bin = 0;     //index of the onset
	double threshold;  // threshold
	
	TH1D* h2 = (TH1D*) histpdf->Clone();
	FilterHistogram(histpdf.GetPtr(),h2, Nfilter);
	
	// Find the maximum after applying the filter
	int maximum = h2->GetMaximum();
	
	// COMPUTE THE THRESHOLD
	//threshold =  Nfilter*(background + ( histpdf->GetMaximum() - background )*fraction);
	threshold =  Nfilter*background + (maximum - Nfilter*background )*fraction;
	
	// nearest and cmp
	// used to find the frequency which is the closest to the threshold, in the particular case no bin is over the threshold
	// this is not used in the analysis
	
	int nearest = 0;
	double cmp = 0;
	
	for(int i = 1; i <= histpdf->GetNbinsX() - Nfilter; i++){ // start the scan
		double sum = 0;
		// apply the filter
		for(int j = 0; j < Nfilter; j++){
			sum += histpdf->GetBinContent(i + j);
		}	
		// find the closes to the threshold
		if(sum > cmp && sum < threshold){
			cmp = sum;
			nearest = static_cast<int>(i);
		}
		// threshold criteria
		if(sum > threshold){
			//std::cout << "identified correct bin " << i << std::endl;
			onset = histpdf->GetBinCenter(static_cast<int>(i));
			bin = i;
			break;
		}
		// uncomment the following two lines to return the closest bin to the threshold, in case no onset is found
		//onset = histpdf->GetBinCenter(nearest); 
		//bin = nearest;
		
		// in case no onset is found, return the last bin as the onset
		onset = histpdf->GetBinCenter(i);
		bin = i;		
	}
	// constan fraction
	std::cout << "Constant fraction:" << std::endl;
	std::cout << "identified frequency: " << onset << " number of bin: " <<  histpdf->GetNbinsX() << std::endl;
	return onset;

}

////////////////////////////////////////
// SIGNIFICANCE

double sumNeighbors(ROOT::RDF::RResultPtr<TH1D> histpdf, double Nsigma, double background, int Nfilter){
	
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	double threshold;
	threshold = Nsigma*Nsigma + 2*Nsigma*sqrt(Nfilter*background);
	for(int i = 1; i < histpdf->GetNbinsX() - 1; i++){
		double sum = 0;
		for(int j = 0; j < Nfilter; j++){
			sum += histpdf->GetBinContent(i + j);
		}
		if(sum > threshold){
			onset = histpdf->GetBinCenter(static_cast<int>(i));
			bin = i;
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	return onset;
}

double Significance(ROOT::RDF::RResultPtr<TH1D> histpdf, double Nsigma, double background, int Nfilter){
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	double threshold;
	
	TH1D* h2 = (TH1D*) histpdf->Clone();
	
	FilterHistogram(histpdf.GetPtr(),h2, Nfilter);
	
	// create the filtered histogram
	std::cout << "check filter: " << std::endl;
	//for(int i = 1; i < histpdf->GetNbinsX() - 1; i++){
	//	std::cout << h2->GetBinContent(i) << std::endl;
	//}

	//threshold =  Nfilter*background + Nsigma*sqrt(Nfilter)*sqrt(background);
	threshold = Nsigma*Nsigma + 2*Nsigma*sqrt(Nfilter*background);
	for(int i = 1; i < histpdf->GetNbinsX() - 1; i++){
		if(h2->GetBinContent(i) > threshold){
			onset = histpdf->GetBinCenter(static_cast<int>(i));
			bin = i;
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	return onset;
}

////////////////////////////////////////
// RUNNING DIFFERENCE
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
	return onset;
}

#endif
