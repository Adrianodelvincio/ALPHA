#include <TMath.h>
#include <cmath>


double linear(double x, double start = 0, double peak = 0,double end = 0){
	//Condizioni raccordo
	double a = (start - peak)/(end - peak);
	double b = (peak - start)/(1 - peak/end);
	if(x <= start){ return 0;}
	else if(x > start && x <= peak){ return (x - start);}
	else if(x > peak && x <= end) { return a*x + b;}
	else { return 0;}
}

double parabola(double x, double start = 0, double peak = 0, double end = 0){
	if(x <= start){ return 0;}
	else if( x > start && x <= peak){ return pow(x - start,2);}
	else if( x > peak && x<= end) {
	double m = -pow(peak - start,2)/(end - peak);
	return m*(x - peak) + pow(peak - start,2);}
	else{return 0;}
}

double Cruijff(double x, double onset, double x0, double sigma0, double sigma1, double k0, double k1, double Norm){
        double result;
        double arg2 = x - x0;
        if(x < onset){
        	result = 0;
        }
        else if(arg2 <= 0 && x >= onset){	// Left tail  
                result = TMath::Exp(-TMath::Power(arg2,2) / (2*sigma0*sigma0 + k0 *TMath::Power(arg2,2)));
        }
        else{					// Right tail
                result = TMath::Exp(-TMath::Power(arg2,2) / (2*sigma1*sigma1 + k1 *TMath::Power(arg2,2)));
        }
        return Norm*result;
}
