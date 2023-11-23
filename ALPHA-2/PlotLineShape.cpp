#include <iostream>
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include <TMath.h>
#include "Headers/toyLineShape.h"

void PlotLineShape(int NBin = 30){
	int Nbin = NBin;
	TNtuple file_pdf1("pdf1", "pdf1","x:y");
	TNtuple file_pdf2("pdf2", "pdf2","x:y");
	file_pdf1.ReadFile("LineShape/lineShape1.csv");
	file_pdf2.ReadFile("LineShape/lineShape2.csv");
	gInterpreter->GenerateDictionary("ToyLine", "Headers/toyLineShape.h");
	vector<double> v1; // Frequency pd1
	vector<double> v2; // Counts per frequence pdf1
	vector<double> t1; // Frequence pdf2
	vector<double> t2; // Counts per frequence pdf2

	for (int row = 0; row < file_pdf1.GetEntries(); ++row) {
		file_pdf1.GetEntry(row);
		v1.push_back(file_pdf1.GetArgs()[0]); 	// Extract frequence
		v2.push_back(file_pdf1.GetArgs()[1]);	// Extract Nls
		}

	for (int row = 0; row < file_pdf2.GetEntries(); ++row) {
		file_pdf2.GetEntry(row);
		t1.push_back(file_pdf2.GetArgs()[0]); 	// Extract frequence
		t2.push_back(file_pdf2.GetArgs()[1]);	// Extract Nls
		}

	Double_t frequence[v1.size()]; Double_t pdf1[v2.size()];
	Double_t frequence2[t1.size()]; Double_t pdf2[t2.size()];

	std::copy(v1.begin(),v1.end(),frequence); 		
	std::copy(v2.begin(),v2.end(),pdf1);
	std::copy(t1.begin(),t1.end(),frequence2); 
	std::copy(t2.begin(),t2.end(),pdf2);

	// Interpolate the data with spline
	TSpline3 *spline1 = new TSpline3("LineShape1", frequence, pdf1, v1.size());
	TSpline3 *spline2 = new TSpline3("LineShape2", frequence2, pdf2, t1.size());
	TH1F *histpdf1 = new TH1F("hist1", "pdf1", Nbin,frequence[0], frequence[v1.size() -1]);
	TH1F *histpdf2 = new TH1F("hist2", "pdf2", Nbin,frequence2[0], frequence2[t1.size() -1]);

	// Set Content Histograms
	SetContent(histpdf1,Nbin,spline1);
	SetContent(histpdf2,Nbin,spline2);
	// Normalize histograms
	Double_t factor = 1.;
	histpdf1->Scale(factor/histpdf1->GetEntries());
	histpdf2->Scale(factor/histpdf2->GetEntries());

	// Fit the data with a function
	auto g1 = new TGraph(v1.size(), frequence,pdf1);
	auto g2 = new TGraph(t1.size(), frequence2,pdf2); 
	// Try convolution of a Gauss pdf and an exponential
	
	TF1 *myExp = new TF1("myExp","expo", 1419, 1425); // build the model
	myExp->SetRange(1419,1425);
	
	TF1 *myGaus = new TF1("myGaus","gaus", 1419, 1425);
	TF1 *f1;
	TF1 *f2;
	TF1 *f3;
	
	//myGaus->SetParameters(13.87,1420.89,0.889);
	//myExp->SetParameters(16,1420.5, 0.8);
	g2->Fit("myExp","","", 1420.5,1425);
	g2->Fit("myGaus");
	
	TCanvas *a = new TCanvas("a", "Pdf 1 model Gauss*Exp");
	auto pad2 = new TPad("pad2", "pad",0,0,1,1);
	pad2->Divide(2,1,0.001,0.001); pad2->Draw();
	pad2->cd(1);
	gStyle->SetOptFit(1);
	g1->SetMarkerStyle(21);
	g1->SetMarkerColor(kBlue);
	g1->SetTitle("Pdf 1");
	g1->GetYaxis()->SetTitle("Counts");
	g1->GetXaxis()->SetTitle("frequencies [MHz]");
	g1->Draw();
	pad2->cd(2);
	g2->SetMarkerStyle(21);
	g2->SetMarkerColor(kBlue);
	g2->SetMarkerStyle(21);
	g2->SetTitle("Pdf 2");
	g2->GetYaxis()->SetTitle("Counts");
	g2->GetXaxis()->SetTitle("frequencies [MHz]");
	g2->Draw();
	myGaus->SetLineColor(kGreen);
	myExp->SetLineColor(kPink);
	myGaus->Draw("same");
	myExp->Draw("same");
	
	/*
	//Visualize spline
	TCanvas *a1 = new TCanvas("a1","Spline");
	auto pad = new TPad("pad", "pad",0,0,1,1);
	pad->Divide(2,1,0.001,0.001); pad->Draw();
	pad->cd(1);
	file_pdf1.SetMarkerStyle(21);
	//file_pdf1.SetTitle("Data Pdf 1");
	file_pdf1.Draw("y:x");
	spline1->SetLineColor(kRed);
	spline1->SetLineWidth(3);
	spline1->Draw("same");
	pad->cd(2);
	//file_pdf2.SetTitle("Data Pdf 2");
	file_pdf2.SetMarkerStyle(21);
	file_pdf2.Draw("y:x");
	spline2->SetLineColor(kRed);
	spline2->SetLineWidth(3);
	spline2->Draw("same");


	TCanvas *a2 = new TCanvas("a2", "Extracted histogram");
	auto pad1 = new TPad("pad1", "pad",0,0,1,1);
	pad1->Divide(2,1,0.001,0.001); pad1->Draw();
	pad1->cd(1);
	histpdf1->Draw();
	pad1->cd(2);
	histpdf2->Draw();*/
}
