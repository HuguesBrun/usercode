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



using namespace RooFit;

void performFits_Bg(){
	// Signal component (Gaussian)
	RooRealVar mgg("mgg","m_{#gamma#gamma}",100.,150., "GeV") ;
	RooArgList observables(mgg); // variables to be generated



	///////////
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
    
    
    cout << "DATACARD float nbBg_cat0 = " << myTreeUnbinned_cat0->GetEntries("mgg>100&&mgg<150") << ";" << endl;
    cout << "DATACARD float nbBg_cat1 = " << myTreeUnbinned_cat1->GetEntries("mgg>100&&mgg<150") << ";" << endl;
    cout << "DATACARD float nbBg_cat2 = " << myTreeUnbinned_cat2->GetEntries("mgg>100&&mgg<150") << ";" << endl;
    cout << "DATACARD float nbBg_cat3 = " << myTreeUnbinned_cat3->GetEntries("mgg>100&&mgg<150") << ";" << endl;


	
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
	return;	


}

	
