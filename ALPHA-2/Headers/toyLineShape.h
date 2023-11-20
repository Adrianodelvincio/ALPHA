#ifndef FITALLFREQUENCY_H
#define FITALLFREQUENCY_H
#include <iostream>
#include <string>
#include "TSystem.h"
#include "TH1.h"
#include "TSpline.h"

void SetContent(TH1 * histpdf, int Nbin, TSpline3 * spline){
	for(int i = 1; i < Nbin; ++i){
		if(spline->Eval(histpdf->GetBinCenter(i)) > 0.){
		histpdf->SetBinContent(i,spline->Eval(histpdf->GetBinCenter(i)));}
		else{
		histpdf->SetBinContent(i,0.);}
		}
}


#endif
