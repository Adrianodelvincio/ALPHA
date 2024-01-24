#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include <ROOT/RDataFrame.hxx>


/////////////////////////////////////////////
// OVER THRESHOLD

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

/////////////////////////////////////////////
// ALGORITHM 2017

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
	//std::cout << "2017: " << " frequence " << histpdf->GetBinCenter(bin)  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

double algorithm_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background, double thr1, double thr2){ //WITH COSMIC BACKGROUND
	// bin 0 is underflow
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
	//std::cout << "2017: " << " i: " << bin  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

/////////////////////////////////////////////
// REVERSED 2017

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

double reverse_2017(ROOT::RDF::RResultPtr<TH1D> histpdf, double background, double Thr1, double Thr2){ // WITH COSMIC BACKGROUND
	// bin 0 is underflow
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
	//std::cout << "2017 reversed" << " i: " << bin << " frequency: " << histpdf->GetBinCenter(bin)  << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	return onset;
}

/////////////////////////////////////////////
// CONSTANT FRACTION ALGORITHM

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

double constFrac(ROOT::RDF::RResultPtr<TH1D> histpdf, double fraction, double background, int Nfilter){ // WITH BACKGROUND
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	double threshold;
	threshold =  Nfilter*(background + ( histpdf->GetMaximum() - background )*fraction);
	
	int nearest;
	double cmp = 0;
	for(int i = 1; i <= histpdf->GetNbinsX() - Nfilter; i++){
		double sum = 0;
		for(int j = 0; j < Nfilter; j++){
			sum += histpdf->GetBinContent(i + j);
			//std::cout << histpdf->GetBinContent(i + j) << " ";
		}
		//std::cout << "		sum: " << sum << std::endl;
		
		if(sum > cmp){
			cmp = sum;
			nearest = static_cast<int>(i + Nfilter/2);
		}
		
		if(sum > threshold){
			if(Nfilter%2 == 0){ 
				onset = histpdf->GetBinCenter(static_cast<int>(i + Nfilter/2));
				bin = static_cast<int>(i + Nfilter/2);
			}
			else{
				onset = histpdf->GetBinCenter(static_cast<int>(i + Nfilter/2 + 0.5));
				bin = static_cast<int>(i + Nfilter/2 + 0.5);
			}
			break;
		}
		
		//onset = histpdf->GetBinCenter(nearest);
		//bin = nearest;
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	//std::cout << "Max: " << histpdf->GetMaximum() << " ConstFrac: " << threshold << " onset: " << onset << " bin content: " << histpdf->GetBinContent(bin) << std::endl;
	//std::cout << "------------" << std::endl;
	return onset;

}

////////////////////////////////////////
//SUMNEIGHBORS OR T-STUDENT TEST

double sumNeighbors(ROOT::RDF::RResultPtr<TH1D> histpdf, double Nsigma, double background, int Nfilter){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	double threshold;
	//threshold =  Nfilter*background + Nsigma*sqrt(Nfilter)*sqrt(background);
	threshold = Nfilter*background + Nsigma*Nsigma + 2*Nsigma*sqrt(Nfilter*background);
	for(int i = 1; i < histpdf->GetNbinsX() - 1; i++){
		double sum = 0;
		for(int j = 0; j < Nfilter; j++){
			sum += histpdf->GetBinContent(i + j);
		}
		if(sum > threshold){
			if(Nfilter%2 == 0){ 
				onset = histpdf->GetBinCenter(static_cast<int>(i + Nfilter/2));
				bin = static_cast<int>(i + Nfilter/2);
			}
			else{
				onset = histpdf->GetBinCenter(static_cast<int>(i + Nfilter/2 + 0.5));
				bin = static_cast<int>(i + Nfilter/2 + 0.5);
			}
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	//std::cout << "sumNeighbors: " << threshold << " background: " << background << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << " bin: " << bin << std::endl;
	return onset;
}

double Significance(ROOT::RDF::RResultPtr<TH1D> histpdf, double Nsigma, double background, int Nfilter){
	// bin 0 is underflow
	double onset = 0;	// onset value 
	double bin = 0;		// bin onset
	double threshold;
	//threshold =  Nfilter*background + Nsigma*sqrt(Nfilter)*sqrt(background);
	threshold = Nsigma*Nsigma + 2*Nsigma*sqrt(Nfilter*background);
	for(int i = 1; i < histpdf->GetNbinsX() - 1; i++){
		double sum = 0;
		for(int j = 0; j < Nfilter; j++){
			sum += histpdf->GetBinContent(i + j);
		}
		if(sum > threshold){
			if(Nfilter%2 == 0){ 
				onset = histpdf->GetBinCenter(static_cast<int>(i + Nfilter/2));
				bin = static_cast<int>(i + Nfilter/2);
			}
			else{
				onset = histpdf->GetBinCenter(static_cast<int>(i + Nfilter/2 + 0.5));
				bin = static_cast<int>(i + Nfilter/2 + 0.5);
			}
			break;
		}
		onset = histpdf->GetBinCenter(i);
		bin = i;
	}
	//std::cout << "sumNeighbors: " << threshold << " background: " << background << " frequency: " << histpdf->GetBinCenter(bin) << " bin content: " << histpdf->GetBinContent(bin) << " bin: " << bin << std::endl;
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
	//std::cout << "runningDiff: " <<  " frequency " << histpdf->GetBinCenter(bin) << std::endl;
	return onset;
	

}

#endif
