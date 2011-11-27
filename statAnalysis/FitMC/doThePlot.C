#include "setTDRStyle.C"


void plotTheHisto(TString summer11, TString spring11){
    
    gStyle->SetOptStat(0);
    
	TFile *myFile = new TFile("huguesOutFile.root");
	TH1D *hSummer = (TH1D*) myFile->Get(summer11);
    hSummer->Sumw2();
	TH1D *hSpring = (TH1D*) myFile->Get(spring11);

    float coeffSummer = 1.0/hSummer->GetEntries();
    float coeffSpring = 1.0/hSpring->GetEntries();
   
    hSummer->Scale(coeffSummer);
    hSpring->Scale(coeffSpring);
    
    TCanvas *c0 = new TCanvas("c0","",600,600);
    c0->SetFillColor(0);
    
    hSpring->SetFillColor(kRed-7w); 
    hSpring->GetXaxis()->SetTitle("R_{9}");
    hSpring->SetTitleOffset(0.8,"X"); 
    hSpring->SetTitleSize(0.05,"X");
    hSpring->SetLabelSize(0.04,"X");
    hSpring->GetYaxis()->SetTitle("nb of #gamma");
	hSpring->SetTitleOffset(1.1,"Y");
    hSpring->SetLabelSize(0.03,"Y");
    hSpring->Draw();    
    hSummer->SetLineColor(kBlue);
   // hSummer->SetLineWidth(3);
    hSummer->Draw("same:E1");
    
    
    TLegend *t = new TLegend(0.16,0.66,0.62,0.86);
    t->SetFillColor(0);
    t->AddEntry(hSpring, "Spring11 MC diphoton born", "F");
    t->AddEntry(hSummer, "Summer11 MC diphoton born", "L");
    t->Draw();
    
    c0->Print(summer11+".gif");
}

doThePlot(){
  	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	/*gROOT->ProcessLine(".x setTDRStyle.C");    */
    
    plotTheHisto("hSummer11_EB","hSpring11_EB");
    plotTheHisto("hSummer11_EE","hSpring11_EE");
    plotTheHisto("hSummer11Loose_EB","hSpring11_EB");
    plotTheHisto("hSummer11Loose_EE","hSpring11_EE");
    plotTheHisto("hSummer11No_EB","hSpring11_EB");
    plotTheHisto("hSummer11No_EE","hSpring11_EE");
    
}
