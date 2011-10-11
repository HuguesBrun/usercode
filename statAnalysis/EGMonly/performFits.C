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

//////////// cut NN
//float nbHiggs_cat0 = 4.95352; // après cut, les 2 dans barrel les > 0.94
//float nbBg_cat0 = 1240.14;

//float nbHiggs_cat1 = 4.78652; // les 2 dans barrel un avec < 0.94
//float nbBg_cat1 = 1824.72;

//float nbHiggs_cat2 = 2.08953; 
//float nbBg_cat2 = 1094.36;

//float nbHiggs_cat3 = 2.30082; 
//float nbBg_cat3 = 1699.16;




/*
 
///////////////////   Cuts optimisés  
 float nbHiggs_cat0 = 5.56407; // après cut, les 2 dans barrel les > 0.94
float nbBg_cat0 = 1509.46;

float nbHiggs_cat1 = 5.4925; // les 2 dans barrel un avec < 0.94
float nbBg_cat1 = 2876.49;

float nbHiggs_cat2 = 2.39515; 
float nbBg_cat2 = 1438.41;

float nbHiggs_cat3 = 2.49775; 
float nbBg_cat3 = 2374.22;*/

/*
///////////////////   EGM 006 

float nbHiggs_cat0 = 3.81186; // après cut, les 2 dans barrel les > 0.94
float nbBg_cat0 = 899.865;

float nbHiggs_cat1 = 4.71195; // les 2 dans barrel un avec < 0.94
float nbBg_cat1 = 2134.1;

float nbHiggs_cat2 = 1.79174; 
float nbBg_cat2 = 945.46;

float nbHiggs_cat3 = 2.30574; 
float nbBg_cat3 = 2198.94;
*/

using namespace RooFit;

void performFits(){
	// Signal component (Gaussian)
	RooRealVar mgg("mgg","m_{#gamma#gamma}",100.,150., "GeV") ;
	RooRealVar r("r","r",1.0,0,10);
	RooArgList observables(mgg); // variables to be generated

	
	TFile *myFile = new TFile("diphoFile_GluGluToHToGGM115.root");

	////// cat 0 /////////
	
	// Signal component (Gaussian)
	TH1F *histoDATA_cat0 = (TH1F*) myFile->Get("higgs2_loose1_cat0_TwoReal");
	TH1F *histoDATA_cat1 = (TH1F*) myFile->Get("higgs2_loose1_cat1_TwoReal");
	TH1F *histoDATA_cat2 = (TH1F*) myFile->Get("higgs2_loose1_cat2_TwoReal");
	TH1F *histoDATA_cat3 = (TH1F*) myFile->Get("higgs2_loose1_cat3_TwoReal");

	cout << histoDATA_cat0->GetEntries() << endl;

	RooRealVar mHiggs_cat0("mHiggs_cat0","m_{0}",116,110,140, "GeV") ;
	RooRealVar wHiggs_cat0("wHiggs_cat0","#sigma_{0}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss_cat0("theGauss_cat0","higgs signal",mgg,mHiggs_cat0,wHiggs_cat0) ;
	RooRealVar nGauss_cat0("nGauss_cat0","N_{0}", 25000,0,50000);
	RooExtendPdf etheGauss_cat0("etheGauss_cat0","higgs signal extended",theGauss_cat0, nGauss_cat0);
	RooRealVar mHiggs2_cat0("mHiggs2_cat0","m_{1}",116,110,140, "GeV") ;
	RooRealVar wHiggs2_cat0("wHiggs2_cat0","#sigma_{1}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss2_cat0("theGauss2_cat0","higgs signal",mgg,mHiggs2_cat0,wHiggs2_cat0) ;
	RooRealVar nGauss2_cat0("nGauss2_cat0","N_{1}", 25000,0,50000);
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
	RooRealVar mHiggs_cat1("mHiggs_cat1","m_{0}",116,110,140, "GeV") ;
	RooRealVar wHiggs_cat1("wHiggs_cat1","#sigma_{0}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss_cat1("theGauss_cat1","higgs signal",mgg,mHiggs_cat1,wHiggs_cat1) ;
	RooRealVar nGauss_cat1("nGauss_cat1","N_{0}", 25000,0,50000);
	RooExtendPdf etheGauss_cat1("etheGauss_cat1","higgs signal extended",theGauss_cat1, nGauss_cat1);
	RooRealVar mHiggs2_cat1("mHiggs2_cat1","m_{1}",116,110,140, "GeV") ;
	RooRealVar wHiggs2_cat1("wHiggs2_cat1","#sigma_{1}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss2_cat1("theGauss2_cat1","higgs signal",mgg,mHiggs2_cat1,wHiggs2_cat1) ;
	RooRealVar nGauss2_cat1("nGauss2_cat1","N_{1}", 25000,0,50000);
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
	RooRealVar mHiggs_cat2("mHiggs_cat2","m_{0}",116,110,140, "GeV") ;
	RooRealVar wHiggs_cat2("wHiggs_cat2","#sigma_{0}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss_cat2("theGauss_cat2","higgs signal",mgg,mHiggs_cat2,wHiggs_cat2) ;
	RooRealVar nGauss_cat2("nGauss_cat2","N_{0}", 25000,0,50000);
	RooExtendPdf etheGauss_cat2("etheGauss_cat2","higgs signal extended",theGauss_cat2, nGauss_cat2);
	RooRealVar mHiggs2_cat2("mHiggs2_cat2","m_{1}",116,110,140, "GeV") ;
	RooRealVar wHiggs2_cat2("wHiggs2_cat2","#sigma_{1}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss2_cat2("theGauss2_cat2","higgs signal",mgg,mHiggs2_cat2,wHiggs2_cat2) ;
	RooRealVar nGauss2_cat2("nGauss2_cat2","N_{1}", 25000,0,50000);
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
	RooRealVar mHiggs_cat3("mHiggs_cat3","m_{0}",116,110,140, "GeV") ;
	RooRealVar wHiggs_cat3("wHiggs_cat3","#sigma_{0}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss_cat3("theGauss_cat3","higgs signal",mgg,mHiggs_cat3,wHiggs_cat3) ;
	RooRealVar nGauss_cat3("nGauss_cat3","N_{0}", 25000,0,50000);
	RooExtendPdf etheGauss_cat3("etheGauss_cat3","higgs signal extended",theGauss_cat3, nGauss_cat3);
	RooRealVar mHiggs2_cat3("mHiggs2_cat3","m_{1}",116,110,140, "GeV") ;
	RooRealVar wHiggs2_cat3("wHiggs2_cat3","#sigma_{1}",1.,0.5,5, "GeV") ;
	RooGaussian theGauss2_cat3("theGauss2_cat3","higgs signal",mgg,mHiggs2_cat3,wHiggs2_cat3) ;
	RooRealVar nGauss2_cat3("nGauss2_cat3","N_{1}", 25000,0,50000);
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
	
	theCanvas->Print("gif/fitSignal.gif");
	
	delete theCanvas;
	delete myFile;
	


	///////////
	TFile *myFile2 = new TFile("diphoFile_DATA.root");
	TFile *myFileUnbinned = new TFile("theUnbinnedDATA.root");
	TTree *myTreeUnbinned_cat0 = (TTree*) myFileUnbinned->Get("higgsCat0");
	RooDataSet *myDataUnbinned_cat0 = new RooDataSet();
	myDataUnbinned_cat0 = new RooDataSet("myDataUnbinned_cat0", "myDataUnbinned_cat0", myTreeUnbinned_cat0, mgg);

	TTree *myTreeUnbinned_cat1 = (TTree*) myFileUnbinned->Get("higgsCat1");
	RooDataSet *myDataUnbinned_cat1 = new RooDataSet();
	myDataUnbinned_cat1 = new RooDataSet("myDataUnbinned_cat1", "myDataUnbinned_cat1", myTreeUnbinned_cat1, mgg);

	TTree *myTreeUnbinned_cat2 = (TTree*) myFileUnbinned->Get("higgsCat2");
	RooDataSet *myDataUnbinned_cat2 = new RooDataSet();
	myDataUnbinned_cat2 = new RooDataSet("myDataUnbinned_cat2", "myDataUnbinned_cat2", myTreeUnbinned_cat2, mgg);

	TTree *myTreeUnbinned_cat3 = (TTree*) myFileUnbinned->Get("higgsCat3");
	RooDataSet *myDataUnbinned_cat3 = new RooDataSet();
	myDataUnbinned_cat3 = new RooDataSet("myDataUnbinned_cat3", "myDataUnbinned_cat3", myTreeUnbinned_cat3, mgg);


	
	////////////////////////////////	
	RooRealVar a_cat0("a_cat0", "a", 10, 0, 50);
	RooRealVar b_cat0("b_cat0", "b", 1, 0, 50);
	RooRealVar c_cat0("c_cat0", "c", 1, 0, 50);
	RooBernstein poly_cat0("poly_cat0", "model2", mgg, RooArgList(a_cat0, b_cat0, c_cat0) ); 
	RooFitResult* r2_cat0 = poly_cat0.fitTo(*myDataUnbinned_cat0, Save());
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float a_cat0 = " << a_cat0.getVal() << ";" << endl;
	cout << "DATACARD float b_cat0 = " << b_cat0.getVal() << ";" << endl;
	cout << "DATACARD float c_cat0 = " << c_cat0.getVal() << ";" << endl;

	RooPlot* frame0 = mgg.frame(Title("cat0"), Bins(50)) ;
	myDataUnbinned_cat0->plotOn(frame0, Name("data"));
	poly_cat0.plotOn(frame0, VisualizeError(*r2_cat0,2, kFALSE), FillColor(kGreen));
	poly_cat0.plotOn(frame0, VisualizeError(*r2_cat0,1, kFALSE), FillColor(kYellow));
	poly_cat0.plotOn(frame0, LineWidth(2), LineColor(kRed), Name("model"));
	myDataUnbinned_cat0->plotOn(frame0);
	
	////////////////////////////////	
	RooRealVar a_cat1("a_cat1", "a", 10, 0, 50);
	RooRealVar b_cat1("b_cat1", "b", 1, 0, 50);
	RooRealVar c_cat1("c_cat1", "c", 1, 0, 50);
	RooBernstein poly_cat1("poly_cat1", "model2", mgg, RooArgList(a_cat1, b_cat1, c_cat1) ); 
	RooFitResult* r2_cat1 = poly_cat1.fitTo(*myDataUnbinned_cat1, Save());
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float a_cat1 = " << a_cat1.getVal() << ";" << endl;
	cout << "DATACARD float b_cat1 = " << b_cat1.getVal() << ";" << endl;
	cout << "DATACARD float c_cat1 = " << c_cat1.getVal() << ";" << endl;
	
	RooPlot* frame1 = mgg.frame(Title("cat1"),Bins(50)) ;
	myDataUnbinned_cat1->plotOn(frame1, Name("data"));
	poly_cat1.plotOn(frame1, VisualizeError(*r2_cat1,2, kFALSE), FillColor(kGreen));
	poly_cat1.plotOn(frame1, VisualizeError(*r2_cat1,1, kFALSE), FillColor(kYellow));
	poly_cat1.plotOn(frame1, LineWidth(2), LineColor(kRed), Name("model"));
	myDataUnbinned_cat1->plotOn(frame1);
	
	////////////////////////////////	
	RooRealVar a_cat2("a_cat2", "a", 10, 0, 50);
	RooRealVar b_cat2("b_cat2", "b", 1, 0, 50);
	RooRealVar c_cat2("c_cat2", "c", 1, 0, 50);
	RooBernstein poly_cat2("poly_cat2", "model2", mgg, RooArgList(a_cat2, b_cat2, c_cat2) ); 
	RooFitResult* r2_cat2 = poly_cat2.fitTo(*myDataUnbinned_cat2, Save());
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float a_cat2 = " << a_cat2.getVal() << ";" << endl;
	cout << "DATACARD float b_cat2 = " << b_cat2.getVal() << ";" << endl;
	cout << "DATACARD float c_cat2 = " << c_cat2.getVal() << ";" << endl;

	RooPlot* frame2 = mgg.frame(Title("cat2"),Bins(50)) ;
	myDataUnbinned_cat2->plotOn(frame2, Name("data"));
	poly_cat2.plotOn(frame2, VisualizeError(*r2_cat2,2, kFALSE), FillColor(kGreen));
	poly_cat2.plotOn(frame2, VisualizeError(*r2_cat2,1, kFALSE), FillColor(kYellow));
	poly_cat2.plotOn(frame2, LineWidth(2), LineColor(kRed), Name("model"));
	myDataUnbinned_cat2->plotOn(frame2);
	
	////////////////////////////////	
	RooRealVar a_cat3("a_cat3", "a", 10, 0, 50);
	RooRealVar b_cat3("b_cat3", "b", 1, 0, 50);
	RooRealVar c_cat3("c_cat3", "c", 1, 0, 50);
	RooBernstein poly_cat3("poly_cat3", "model2", mgg, RooArgList(a_cat3, b_cat3, c_cat3) ); 
	RooFitResult* r2_cat3 = poly_cat3.fitTo(*myDataUnbinned_cat3, Save());
	cout << "DATACARD //////////////////////////////////" << endl;
	cout << "DATACARD float a_cat3 = " << a_cat3.getVal() << ";" << endl;
	cout << "DATACARD float b_cat3 = " << b_cat3.getVal() << ";" << endl;
	cout << "DATACARD float c_cat3 = " << c_cat3.getVal() << ";" << endl;

	RooPlot* frame3 = mgg.frame(Title("cat3"),Bins(50)) ;
	myDataUnbinned_cat3->plotOn(frame3, Name("data"));
	poly_cat3.plotOn(frame3, VisualizeError(*r2_cat3,2, kFALSE), FillColor(kGreen));
	poly_cat3.plotOn(frame3, VisualizeError(*r2_cat3,1, kFALSE), FillColor(kYellow));
	poly_cat3.plotOn(frame3, LineWidth(2), LineColor(kRed), Name("model"));
	myDataUnbinned_cat3->plotOn(frame3);

	/////////////////////////////////////////////
	
	TCanvas *theCanvas = new TCanvas("theCanvas","coucou",1200,1200);
	theCanvas->Divide(2,2);
	theCanvas->cd(1);
	frame0->Draw();	
	RooArgSet *r2_cat0_param = poly_cat0.getVariables();
	plotParameters(r2_cat0_param, theCanvas, 1, frame0, false);
	theCanvas->cd(2);
	frame1->Draw();	
	RooArgSet *r2_cat1_param = poly_cat1.getVariables();
	plotParameters(r2_cat1_param, theCanvas, 2, frame1, false);
	theCanvas->cd(3);
	frame2->Draw();	
	RooArgSet *r2_cat2_param = poly_cat2.getVariables();
	plotParameters(r2_cat2_param, theCanvas, 3, frame2, false);
	theCanvas->cd(4);
	frame3->Draw();	
	RooArgSet *r2_cat3_param = poly_cat3.getVariables();
	plotParameters(r2_cat3_param, theCanvas, 4, frame3, false);
	theCanvas->Print("gif/fitBg.gif");
	
	delete theCanvas;
	delete myFile2;
	return;	


}

	
