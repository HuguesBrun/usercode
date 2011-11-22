#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TLatex.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include "plotParameters.cc"

#pragma optimize 0

TString FloatToString(float number){
	ostringstream oss;
	oss << number;
	return oss.str();
}
float intLumi = 1143;
TString massPoint = "110";
float massPointFloat = 110.;

using namespace RooFit;


TFile *myFileHiggsGG = new TFile("diphoFile_GluGluToHToGGM"+massPoint+".root");
TFile *myFileHiggsWHZHH = new TFile("diphoFile_WHZHHToGGM"+massPoint+".root");
TFile *myFileHiggsTTHH = new TFile("diphoFile_TTHHToGGM"+massPoint+".root");
TFile *myFileHiggsVBFH = new TFile("diphoFile_VBFHToGGM"+massPoint+".root");	



TH1F *histoWithWeight(TString nom){
	TH1F *theHistoGG = (TH1F*) myFileHiggsGG->Get(nom);
	TH1F *theHistoWHZHH = (TH1F*) myFileHiggsWHZHH->Get(nom);
	TH1F *theHistoTTHH= (TH1F*) myFileHiggsTTHH->Get(nom);
	TH1F *theHistoVBFH = (TH1F*) myFileHiggsVBFH->Get(nom);
	
	float weights[4] = {
		16.63*0.00159/109993, // Higgs GLuon fusion
		(0.8754+0.4721)*0.00159/109964, // Higgs strallung
		0.1257*0.00159/22000, // Higgs TT
		1.391*0.00159/109692 // Higgs VBF
    };  
	theHistoGG->Sumw2();
	theHistoGG->Scale(weights[0]);
	theHistoGG->Add(theHistoWHZHH, weights[1]);
	theHistoGG->Add(theHistoTTHH, weights[2]);
	theHistoGG->Add(theHistoVBFH, weights[3]);
	theHistoGG->Scale(intLumi);
	return theHistoGG;
}

void performFits_110(){
	// Signal component (Gaussian)
	RooRealVar mgg("mgg","m_{#gamma#gamma}",100.,150., "GeV") ;
	RooRealVar r("r","r",1.0,0,10);
	RooArgList observables(mgg); // variables to be generated

	


	////// cat 0 /////////
	
	// Signal component (Gaussian)
	TH1F *histoDATA_cat0 = histoWithWeight("higgs2_CICB_cat0_TwoReal");
	TH1F *histoDATA_cat1 = histoWithWeight("higgs2_CICB_cat1_TwoReal");
	TH1F *histoDATA_cat2 = histoWithWeight("higgs2_CICB_cat2_TwoReal");
	TH1F *histoDATA_cat3 = histoWithWeight("higgs2_CICB_cat3_TwoReal");
	
	

	cout << histoDATA_cat0->GetEntries() << endl;

	
	cout << "DATACARD float nbHiggs_cat0 = " << histoDATA_cat0->Integral() << ";" << endl;
	cout << "DATACARD float nbHiggs_cat1 = " << histoDATA_cat1->Integral() << ";" <<endl;
	cout << "DATACARD float nbHiggs_cat2 = " << histoDATA_cat2->Integral() << ";" <<endl;
	cout << "DATACARD float nbHiggs_cat3 = " << histoDATA_cat3->Integral() << ";" <<endl;

	
	
	
	RooRealVar mHiggs_cat0("mHiggs_cat0","m_{0}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs_cat0("wHiggs_cat0","#sigma_{0}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss_cat0("theGauss_cat0","higgs signal",mgg,mHiggs_cat0,wHiggs_cat0) ;
	RooRealVar nGauss_cat0("nGauss_cat0","N_{0}", 25,0,50);
	RooExtendPdf etheGauss_cat0("etheGauss_cat0","higgs signal extended",theGauss_cat0, nGauss_cat0);
	RooRealVar mHiggs2_cat0("mHiggs2_cat0","m_{1}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs2_cat0("wHiggs2_cat0","#sigma_{1}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss2_cat0("theGauss2_cat0","higgs signal",mgg,mHiggs2_cat0,wHiggs2_cat0) ;
	RooRealVar nGauss2_cat0("nGauss2_cat0","N_{1}", 25,0,50);
	RooExtendPdf etheGauss2_cat0("etheGauss2_cat0","higgs signal extended",theGauss2_cat0, nGauss2_cat0);
	
	RooAddPdf model_cat0("model_cat0","model",RooArgList(etheGauss_cat0, etheGauss2_cat0));

	RooDataHist theData_cat0("theData_cat0","data",mgg,histoDATA_cat0);
	model_cat0.fitTo(theData_cat0);
	
	
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float mHiggs_cat0Val = " << mHiggs_cat0.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs_cat0Val = " << wHiggs_cat0.getVal() << ";" << endl;
	cout << "DATACARD float nGauss_cat0 = " << nGauss_cat0.getVal() << ";" << endl;
	cout << "DATACARD float mHiggs2_cat0Val = " << mHiggs2_cat0.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs2_cat0Val = " << wHiggs2_cat0.getVal() << ";" << endl;
	cout << "DATACARD float nGauss2_cat0 = " << nGauss2_cat0.getVal() << ";" << endl;
	
	
	RooPlot* frame1 = mgg.frame(Title("cat0"),Bins(70)) ;
	theData_cat0.plotOn(frame1, Name("mc"));
	model_cat0.plotOn(frame1, Components(etheGauss_cat0),LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
	model_cat0.plotOn(frame1, Components(etheGauss2_cat0),LineStyle(kDashed), LineColor(kGreen), LineWidth(2));
	model_cat0.plotOn(frame1, LineColor(kBlue), LineWidth(2), Name("model"));

	
	//////////// cat 1 
	RooRealVar mHiggs_cat1("mHiggs_cat1","m_{0}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs_cat1("wHiggs_cat1","#sigma_{0}",2.,0.5,6., "GeV") ;
	RooGaussian theGauss_cat1("theGauss_cat1","higgs signal",mgg,mHiggs_cat1,wHiggs_cat1) ;
	RooRealVar nGauss_cat1("nGauss_cat1","N_{0}", 25,0,50);
	RooExtendPdf etheGauss_cat1("etheGauss_cat1","higgs signal extended",theGauss_cat1, nGauss_cat1);
	RooRealVar mHiggs2_cat1("mHiggs2_cat1","m_{1}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs2_cat1("wHiggs2_cat1","#sigma_{1}",2.,0.5,6., "GeV") ;
	RooGaussian theGauss2_cat1("theGauss2_cat1","higgs signal",mgg,mHiggs2_cat1,wHiggs2_cat1) ;
	RooRealVar nGauss2_cat1("nGauss2_cat1","N_{1}", 25,0,50);
	RooExtendPdf etheGauss2_cat1("etheGauss2_cat1","higgs signal extended",theGauss2_cat1, nGauss2_cat1);
	
	RooAddPdf model_cat1("model_cat1","model",RooArgList(etheGauss_cat1, etheGauss2_cat1));
	
	RooDataHist theData_cat1("theData_cat1","data",mgg,histoDATA_cat1);
	model_cat1.fitTo(theData_cat1);
	
	
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float mHiggs_cat1Val = " << mHiggs_cat1.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs_cat1Val = " << wHiggs_cat1.getVal() << ";" << endl;
	cout << "DATACARD float nGauss_cat1 = " << nGauss_cat1.getVal() << ";" << endl;
	cout << "DATACARD float mHiggs2_cat1Val = " << mHiggs2_cat1.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs2_cat1Val = " << wHiggs2_cat1.getVal() << ";" << endl;
	cout << "DATACARD float nGauss2_cat1 = " << nGauss2_cat1.getVal() << ";" << endl;
	
	RooPlot* frame2 = mgg.frame(Title("cat1"),Bins(70)) ;
	theData_cat1.plotOn(frame2, Name("mc"));
	model_cat1.plotOn(frame2, Components(etheGauss_cat1), LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
	model_cat1.plotOn(frame2, Components(etheGauss2_cat1),LineStyle(kDashed), LineColor(kGreen), LineWidth(2));
	model_cat1.plotOn(frame2, LineColor(kBlue), LineWidth(2), Name("model"));
	
	
	//////////// cat 2
	RooRealVar mHiggs_cat2("mHiggs_cat2","m_{0}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs_cat2("wHiggs_cat2","#sigma_{0}",1.,0.3,5, "GeV") ;
	RooGaussian theGauss_cat2("theGauss_cat2","higgs signal",mgg,mHiggs_cat2,wHiggs_cat2) ;
	RooRealVar nGauss_cat2("nGauss_cat2","N_{0}", 25,0,50);
	RooExtendPdf etheGauss_cat2("etheGauss_cat2","higgs signal extended",theGauss_cat2, nGauss_cat2);
	RooRealVar mHiggs2_cat2("mHiggs2_cat2","m_{1}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs2_cat2("wHiggs2_cat2","#sigma_{1}",1.,0.3,5, "GeV") ;
	RooGaussian theGauss2_cat2("theGauss2_cat2","higgs signal",mgg,mHiggs2_cat2,wHiggs2_cat2) ;
	RooRealVar nGauss2_cat2("nGauss2_cat2","N_{1}", 25,0,50);
	RooExtendPdf etheGauss2_cat2("etheGauss2_cat2","higgs signal extended",theGauss2_cat2, nGauss2_cat2);

	RooAddPdf model_cat2("model_cat2","model",RooArgList(etheGauss_cat2, etheGauss2_cat2));
	
	RooDataHist theData_cat2("theData_cat2","data",mgg,histoDATA_cat2);
	model_cat2.fitTo(theData_cat2);
	
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float mHiggs_cat2Val = " << mHiggs_cat2.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs_cat2Val = " << wHiggs_cat2.getVal() << ";" << endl;
	cout << "DATACARD float nGauss_cat2 = " << nGauss_cat2.getVal() << ";" << endl;
	cout << "DATACARD float mHiggs2_cat2Val = " << mHiggs2_cat2.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs2_cat2Val = " << wHiggs2_cat2.getVal() << ";" << endl;
	cout << "DATACARD float nGauss2_cat2 = " << nGauss2_cat2.getVal() << ";" << endl;
	
	RooPlot* frame3 = mgg.frame(Title("cat2"),Bins(70)) ;
	theData_cat2.plotOn(frame3, Name("mc"));
	model_cat2.plotOn(frame3, Components(etheGauss_cat2), LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
	model_cat2.plotOn(frame3, Components(etheGauss2_cat2),LineStyle(kDashed), LineColor(kGreen), LineWidth(2));
	model_cat2.plotOn(frame3, LineColor(kBlue), LineWidth(2), Name("model"));
	
	//////////// cat 3
	RooRealVar mHiggs_cat3("mHiggs_cat3","m_{0}",massPointFloat,massPointFloat-2.,massPointFloat+2.0, "GeV") ;
	RooRealVar wHiggs_cat3("wHiggs_cat3","#sigma_{0}",1.,0.3,5, "GeV") ;
	RooGaussian theGauss_cat3("theGauss_cat3","higgs signal",mgg,mHiggs_cat3,wHiggs_cat3) ;
	RooRealVar nGauss_cat3("nGauss_cat3","N_{0}", 25,0,50);
	RooExtendPdf etheGauss_cat3("etheGauss_cat3","higgs signal extended",theGauss_cat3, nGauss_cat3);
	RooRealVar mHiggs2_cat3("mHiggs2_cat3","m_{1}",massPointFloat,massPointFloat-2.,massPointFloat+2., "GeV") ;
	RooRealVar wHiggs2_cat3("wHiggs2_cat3","#sigma_{1}",1.,0.3,5, "GeV") ;
	RooGaussian theGauss2_cat3("theGauss2_cat3","higgs signal",mgg,mHiggs2_cat3,wHiggs2_cat3) ;
	RooRealVar nGauss2_cat3("nGauss2_cat3","N_{1}", 25,0,50);
	RooExtendPdf etheGauss2_cat3("etheGauss2_cat3","higgs signal extended",theGauss2_cat3, nGauss2_cat3);

	RooAddPdf model_cat3("model_cat3","model",RooArgList(etheGauss_cat3, etheGauss2_cat3));
	
	RooDataHist theData_cat3("theData_cat3","data",mgg,histoDATA_cat3);
	model_cat3.fitTo(theData_cat3);
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float mHiggs_cat3Val = " << mHiggs_cat3.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs_cat3Val = " << wHiggs_cat3.getVal() << ";" << endl;
	cout << "DATACARD float nGauss_cat3 = " << nGauss_cat3.getVal() << ";" << endl;
	cout << "DATACARD float mHiggs2_cat3Val = " << mHiggs2_cat3.getVal() << ";" << endl;
	cout << "DATACARD float wHiggs2_cat3Val = " << wHiggs2_cat3.getVal() << ";" << endl;
	cout << "DATACARD float nGauss2_cat3 = " << nGauss2_cat3.getVal() << ";" << endl;
	
	RooPlot* frame4 = mgg.frame(Title("cat3"),Bins(70)) ;
	theData_cat3.plotOn(frame4, Name("mc"));
	model_cat3.plotOn(frame4, Components(etheGauss_cat3), LineStyle(kDashed), LineColor(kMagenta), LineWidth(2));
	model_cat3.plotOn(frame4, Components(etheGauss2_cat3),LineStyle(kDashed), LineColor(kGreen), LineWidth(2));
	model_cat3.plotOn(frame4, LineWidth(2), Name("model"));
	
	
	
	
	
	
	
	TCanvas *theCanvas = new TCanvas("theCanvas","coucou",1200,1200);
	theCanvas->Divide(2,2);
	theCanvas->cd(1);
	frame1->Draw();
  RooArgSet *signal_cat0_param = model_cat0.getVariables();
  plotParameters(signal_cat0_param, theCanvas, 1, frame1, true);
	theCanvas->cd(2);
	frame2->Draw();
  RooArgSet *signal_cat1_param = model_cat1.getVariables();
  plotParameters(signal_cat1_param, theCanvas, 2, frame2, true);
	theCanvas->cd(3);
	frame3->Draw();
  RooArgSet *signal_cat2_param = model_cat2.getVariables();
  plotParameters(signal_cat2_param, theCanvas, 3, frame3, true);
	theCanvas->cd(4);
	frame4->Draw();
  RooArgSet *signal_cat3_param = model_cat3.getVariables();
  plotParameters(signal_cat3_param, theCanvas, 4, frame4, true);
	
	theCanvas->Print("gif/fitSignal_M"+massPoint+".gif");
	
	delete theCanvas;
	



}

	
