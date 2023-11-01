#include <iostream>
#include <TMath.h>

double Rayleigh(double *, double *);
double TemplateWall(double *, double *);
double TemplateBK(double *, double *);
double ModelFit(double *, double *);

TCanvas *canvas1;

void FitDistribution(){

    // FUNZIONI DI FIT
	TF1 *func1 = new TF1("func1", Rayleigh,0,4,2);
	TF1 *func2 = new TF1("func2", TemplateWall,0,4,1);
	TF1 *func3 = new TF1("func3", TemplateBK,0,4,1);
    TF1 *modelloFit = new TF1("model", ModelFit , 0,4,4);
    Double_t paramsRay[2] = {1,1};
    Double_t paramWall[1] = {1};
    func1->SetParameters(paramsRay);
    func2->SetParameters(paramWall);
    func3->SetParameters(paramWall);
    
    
    // load the data, compute R and cut1
    auto tdf = ROOT::RDF::MakeCsvDataFrame("Dataset/r68465_uw_exp_freq4.vertex.csv");
    auto tdf1 = tdf.Define("Radius", "TMath::Sqrt(X*X + Y*Y)");
    auto tdf2 = tdf1.Display({"CutsType1", "X"}, 10);
    tdf2->Print();
    
    auto Histof4 = tdf1.Filter("CutsType1 ==  \" 1\"").Histo1D({"Counts Frequency 4","Counts",30u,0.,4.}, "Radius");
    
    // Now the fit
    Double_t paramsModel[4] = { 0.5, // Normalizzazione Parete
                           0.5, // Normalizzazione centro
                           0.5, // Sigma Rayleigh
                           1,}; // Cosmici
    modelloFit->SetParameters(paramsModel);
    Histof4->Fit("model", "", "same",0,4);
    TF1 *fitUp = Histof4->GetFunction("model");
    double p1 = fitUp->GetParameter(0);
    double p2 = fitUp->GetParameter(1);
    double p3 = fitUp->GetParameter(2);
    double p4 = fitUp->GetParameter(3);
   
    // Show the Data and the Fit 
    canvas1 = new TCanvas("a", "Rayleigh Distribution", 1100, 900);
    gStyle->SetOptFit(1011);
    //modelloFit->Draw();
    Histof4->DrawClone();
    
    double_t paramsRay1[2] = {p3,p4};
    func1->SetParameters(paramsRay1);
    func1->SetLineColor(kBlue);
	//func1->Draw("same");
    
    double_t paramsWall1[1] = {p1};
    func2->SetParameters(paramsWall1);
    func2->SetLineColor(kGreen);
	//func2->Draw("same");
    
    double_t paramsWall2[1] = {p2};
    func3->SetParameters(paramsWall2);
    func3->SetLineColor(kOrange);
	//func3->Draw("same");
	
    
    canvas1->Modified();
	canvas1->Update();

}

double Rayleigh(double *x, double *par){
        double xx = x[0];
        double Nray = par[1];
        double sigma = par[0];
        return Nray*(xx/(TMath::Power(sigma,2))) * TMath::Exp( -xx*xx/(2*TMath::Power(sigma,2)));}

double TemplateWall(double *x, double *par){
        double xx = x[0];
        double Nwall = par[0];
        //double bins[] = {1,2,3,4,5};
        //double Counts[] = {1,2};
    
        double bins[30] = {0.133, 0.267, 0.4, 0.533, 0.667, 0.8, 0.933,
            1.067, 1.2, 1.333, 1.467, 1.6, 1.733, 1.867, 2.0, 2.133, 2.267,
            2.4, 2.533, 2.667, 2.8, 2.933, 3.067, 3.2, 3.333, 3.467, 3.6,
            3.733, 3.867, 4.0};

        double Counts[30] = {0.007, 0.027, 0.047, 0.064, 0.082, 0.101, 0.126,
        0.136, 0.185, 0.216, 0.256, 0.322, 0.348, 0.409, 0.397, 0.458, 0.479,
        0.484, 0.504, 0.432, 0.382, 0.376, 0.353, 0.297, 0.267, 0.24, 0.206,
        0.114, 0.088, 0.097};

        const int arraysize = 30;

    int index = 0;
        for(int i = 0; i < arraysize; i++){
            if(xx >= bins[i]){
            index = i;
            }
        }

    return Nwall*Counts[index];
}

double TemplateBK(double *x, double *par){
    double xx = x[0];
    double Nbk = par[0];
        double bins[30] = {0.133, 0.267, 0.4, 0.533, 0.667, 0.8, 0.933,
        1.067, 1.2, 1.333, 1.467, 1.6, 1.733, 1.867, 2.0, 2.133, 2.267,
        2.4, 2.533, 2.667, 2.8, 2.933, 3.067, 3.2, 3.333, 3.467, 3.6,
        3.733, 3.867, 4.0};

        double Counts[30] = {0.017, 0.025, 0.021, 0.059, 0.038, 0.097,
        0.08, 0.109, 0.189, 0.155, 0.168, 0.185, 0.185, 0.256, 0.181,
        0.281, 0.269, 0.328, 0.302, 0.344, 0.378, 0.428, 0.433, 0.37,
        0.361, 0.475, 0.588, 0.395, 0.403, 0.382};

    	int arraysize = 30;

    	int index = 0;
    	for(int i = 0; i < arraysize; i++){
        	if(xx >= bins[i]){
            	index = i;
        	}
    	}

    return Nbk*Counts[index];}

double ModelFit(double *x, double *par){
    double xx = x[0];
    double Nwall = par[0];
    double Nray = par[1];
    double sigma = par[3];
    double Nbk = par[4];
    
    /* TEMPLATE WALL*/
    double bins[30] = {0.133, 0.267, 0.4, 0.533, 0.667, 0.8, 0.933,
            1.067, 1.2, 1.333, 1.467, 1.6, 1.733, 1.867, 2.0, 2.133, 2.267,
            2.4, 2.533, 2.667, 2.8, 2.933, 3.067, 3.2, 3.333, 3.467, 3.6,
            3.733, 3.867, 4.0};

        double WallCounts[30] = {0.007, 0.027, 0.047, 0.064, 0.082, 0.101, 0.126,
        0.136, 0.185, 0.216, 0.256, 0.322, 0.348, 0.409, 0.397, 0.458, 0.479,
        0.484, 0.504, 0.432, 0.382, 0.376, 0.353, 0.297, 0.267, 0.24, 0.206,
        0.114, 0.088, 0.097};

        const int arraysize = 30;

    int index = 0;
        for(int i = 0; i < arraysize; i++){
            if(xx >= bins[i]){
            index = i;
            }
        }
    double tmp = WallCounts[index];
    
    /* TEMPLATE RAY*/
    double ray = (xx/(TMath::Power(sigma,2))) * TMath::Exp( -xx*xx/(2*TMath::Power(sigma,2)));
    
    /*TEMPLATE BACKGROUND*/
    double BKCounts[30] = {0.017, 0.025, 0.021, 0.059, 0.038, 0.097,
        0.08, 0.109, 0.189, 0.155, 0.168, 0.185, 0.185, 0.256, 0.181,
        0.281, 0.269, 0.328, 0.302, 0.344, 0.378, 0.428, 0.433, 0.37,
        0.361, 0.475, 0.588, 0.395, 0.403, 0.382};

    	index = 0;
    	for(int i = 0; i < arraysize; i++){
        	if(xx >= bins[i]){
            	index = i;
        	}
    	}
    double cmt = BKCounts[index];
    
    /*RETURN*/
    return abs(Nwall)*tmp + abs(Nray)*ray + abs(Nbk) * cmt;
}
