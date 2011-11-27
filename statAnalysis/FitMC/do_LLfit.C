#include "setTDRStyle.C"

TFile *myFile = new TFile("huguesOutFile.root");



performTheFit(TString nomHisto){
    TH1D *hLL = (TH1D*) myFile->Get(nomHisto);
    
    
    TCanvas *c0 = new TCanvas("c0","coucou", 600, 600);
    c0->SetFillColor(0);
    hLL->GetXaxis()->SetTitle("k = R9_{Spring11}/R9_{Summer11}");
    hLL->SetTitleOffset(1.2,"X"); 
    hLL->SetTitleSize(0.05,"X");
    hLL->SetLabelSize(0.04,"X");
    hLL->GetYaxis()->SetTitle("LL value");
	hLL->SetTitleOffset(1.5,"Y");
    hLL->SetLabelSize(0.03,"Y");
    hLL->SetMarkerSize(0.5);
	hLL->SetMarkerStyle(20);
    hLL->Draw("E1");

    TF1 *f1 = new TF1("f1", "pol2", 0.998,1.01);
    f1->SetLineColor(kBlue);
    f1->SetLineWidth(3);
    hLL->Fit("f1","R");
    hLL->Draw("E1");
    
    c0->Print(nomHisto+".gif");
    
    float maxi = -1.0 * f1->GetParameter(1)/(2*f1->GetParameter(2));
    float error = 1.0/sqrt(-2*f1->GetParameter(2));
    
    cout << "kmax =  " << maxi << " +- " << error << endl;

}



do_LLfit(){
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
    gStyle->SetOptFit(0);
	gROOT->ProcessLine(".x setTDRStyle.C");    
    
    performTheFit("LLbarrel");
    performTheFit("LLendcap");

    
}