#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooProduct.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooSimWSTool.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"

////// numberOfCat
int nbOfCats = 4;

float nbBg[4] = {1240.14, 1824.72, 1094.36, 1699.16};

TString WSsuffixe = "EGMLOOSE006";



/*//////////// cut NN
float nbHiggs_cat0 = 4.95352; // après cut, les 2 dans barrel les > 0.94
float nbBg_cat0 = 1240.14;

float nbHiggs_cat1 = 4.78652; // les 2 dans barrel un avec < 0.94
float nbBg_cat1 = 1824.72;

float nbHiggs_cat2 = 2.08953; 
float nbBg_cat2 = 1094.36;

float nbHiggs_cat3 = 2.30082; 
float nbBg_cat3 = 1699.16;*/

/*///////////////////   Cuts optimisés  
float nbHiggs_cat0 = 5.56407; // après cut, les 2 dans barrel les > 0.94
float nbBg_cat0 = 1509.46;

float nbHiggs_cat1 = 5.4925; // les 2 dans barrel un avec < 0.94
float nbBg_cat1 = 2876.49;

float nbHiggs_cat2 = 2.39515; 
float nbBg_cat2 = 1438.41;

float nbHiggs_cat3 = 2.49775; 
float nbBg_cat3 = 2374.22;*/


/// EGM loose 006 
float nbHiggs_cat0 = 5.00843; // après cut, les 2 dans barrel les > 0.94
float nbBg_cat0 = 837;

float nbHiggs_cat1 = 6.01467; // les 2 dans barrel un avec < 0.94
float nbBg_cat1 = 2164;

float nbHiggs_cat2 = 2.16363; 
float nbBg_cat2 = 861;

float nbHiggs_cat3 = 2.87662; 
float nbBg_cat3 = 2542;

//////////////////////////////////
float mHiggs_cat0Val = 114.45;
float wHiggs_cat0Val = 1.09446;
float nGauss_cat0 = 8722.75;
float mHiggs2_cat0Val = 112.835;
float wHiggs2_cat0Val = 3.40626;
float nGauss2_cat0 = 2316.34;
//////////////////////////////////
float mHiggs_cat1Val = 115.525;
float wHiggs_cat1Val = 1.59937;
float nGauss_cat1 = 10479;
float mHiggs2_cat1Val = 113.857;
float wHiggs2_cat1Val = 4.48608;
float nGauss2_cat1 = 2904.01;
//////////////////////////////////
float mHiggs_cat2Val = 114.533;
float wHiggs_cat2Val = 1.77269;
float nGauss_cat2 = 3438.39;
float mHiggs2_cat2Val = 114.036;
float wHiggs2_cat2Val = 4.14505;
float nGauss2_cat2 = 1322.49;
//////////////////////////////////
float mHiggs_cat3Val = 115.811;
float wHiggs_cat3Val = 2.40102;
float nGauss_cat3 = 4595.33;
float mHiggs2_cat3Val = 114.324;
float wHiggs2_cat3Val = 5;
float nGauss2_cat3 = 1796.6;
//////////////////////////////////
float tau_cat0Val = -0.0297933;
float nbkg_cat0Val = 417.014;
float tau2_cat0Val = -0.0298167;
float nbkg2_cat0Val = 419.984;
//////////////////////////////////
float tau_cat1Val = -0.033053;
float nbkg_cat1Val = 1060.39;
float tau2_cat1Val = -0.0330616;
float nbkg2_cat1Val = 1103.59;
//////////////////////////////////
float tau_cat2Val = -0.0157393;
float nbkg_cat2Val = 430.694;
float tau2_cat2Val = -0.0157312;
float nbkg2_cat2Val = 430.241;
//////////////////////////////////
float tau_cat3Val = -0.013956;
float nbkg_cat3Val = 1264.13;
float tau2_cat3Val = -0.0361855;
float nbkg2_cat3Val = 1277.91;


TChain *chain0 = new TChain("higgsCat0");
TChain *chain1 = new TChain("higgsCat1");
TChain *chain2 = new TChain("higgsCat2");
TChain *chain3 = new TChain("higgsCat3");


using namespace RooFit;
//using namespace RooStats;

void prepWorkSpace_M115_NN(){
	RooRandom::randomGenerator()->SetSeed(3007);
	
	
	chain0->Add("theUnbinnedDATA.root");
	chain1->Add("theUnbinnedDATA.root");
	chain2->Add("theUnbinnedDATA.root");
	chain3->Add("theUnbinnedDATA.root");
	
/*	TH1F *histoDATAbg_cat0 = (TH1F*) myFile->Get("higgs2_loose1_cat0");
	TH1F *histoDATAbg_cat1 = (TH1F*) myFile->Get("higgs2_loose1_cat1");
	TH1F *histoDATAbg_cat2 = (TH1F*) myFile->Get("higgs2_loose1_cat2");
	TH1F *histoDATAbg_cat3 = (TH1F*) myFile->Get("higgs2_loose1_cat3");*/
	
	

	RooRealVar mgg("mgg","",100.,150.) ;
	RooRealVar r("r","r",1.0,0.,10.);
	//RooArgList observables(mgg); // variables to be generated

	
	///////////////////////////////////////////////////////////////////////////	
	// Signal component (Gaussian)

	RooRealVar mHiggs("mHiggs","m(higgs)",115.,110.,150.) ;
	RooRealVar wHiggs("wHiggs","m(higgs) resolution",0.1,0.,5. ) ;
	RooGaussian theGauss("theGauss","higgs signal",mgg,mHiggs,wHiggs) ;
	
	
	RooRealVar mHiggs2("mHiggs2","m(higgs)",115.,110.,150.) ;
	RooRealVar wHiggs2("wHiggs2","m(higgs) resolution",0.1,0.,5. ) ;
	RooGaussian theGauss2("theGauss2","higgs signal",mgg,mHiggs2,wHiggs2) ;
	
	RooRealVar nbGauss("nbGauss","nb gauss", 10000.,0.,100000.);
	RooRealVar nbGauss2("nbGauss2","nb gauss", 10000.,0.,100000.);
	RooAddPdf sig_pdf("sig_pdf","",RooArgList(theGauss,theGauss2),RooArgList(nbGauss,nbGauss2));
	
	
	
	// bg in double exponential 
	RooRealVar tau("tau","coefficient Tau",-0.1,-2.,0.);
	RooExponential theExp("theExp","background pdf",mgg,tau);
	RooRealVar tau2("tau2","coefficient Tau",-0.1,-2.,0.);
	RooExponential theExp2("theExp2","background pdf",mgg,tau2);
	
	
	RooRealVar nbExp("nbExp","proportion", 300.,0.,10000.);
	RooRealVar nbExp2("nbExp2","proportion", 300.,0.,10000.);
	RooAddPdf bkg_pdf("bkg_pdf","",RooArgList(theExp,theExp2),RooArgList(nbExp,nbExp2));
	
	/// now Sig + Bg model 
	RooRealVar theNbHiggs("theNbHiggs","theNbHiggs",10.,0.,100.);
	RooProduct sig_yield("sig_yield","sig_yield", RooArgSet(theNbHiggs,r));
	RooRealVar bkg_yield("bkg_yield","bkg_yield",300.,0.,30000.);
	RooExtendPdf bkg_ext_pdf("bkg_ext_pdf","bkg_ext_pdf",bkg_pdf,bkg_yield);	// here the model of bg only 
	RooAddPdf tot_pdf("tot_pdf","",RooArgList(sig_pdf,bkg_pdf),RooArgList(sig_yield,bkg_yield));
	
	
	/*RooRealVar bkg_yield_syst("bkg_yield_syst","bkg_yield_syst",300,0,30000);
	RooRealVar bkg_yield_width("bkg_yield_width","bkg_yield_width",12,0,30000);
	RooGaussian bkg_yield_prior("bkg_yield_prior","",bkg_yield,bkg_yield_syst,bkg_yield_width);*/

	
	// now define the categories 
	
	RooCategory cats4("cats4","cats4") ;
	cats4.defineType("cat0") ;
	cats4.defineType("cat1") ;	
	cats4.defineType("cat2") ;
	cats4.defineType("cat3") ;
	
	RooArgList observables(mgg,cats4); // variables to be generated
	
	RooWorkspace myWS("myWS");
	myWS.import(RooArgSet(tot_pdf,cats4)) ;
	//myWS.import(RooArgSet(bkg_yield_prior,cats4));
	myWS.Print("v");
	
	RooSimWSTool sct(myWS) ;
	
	RooSimultaneous* tot_pdf_sim = sct.build("tot_pdf_sim","tot_pdf",SplitParam("mHiggs,wHiggs,nbGauss,mHiggs2,wHiggs2,nbGauss2,tau,nbExp,tau2,nbExp2,bkg_yield,theNbHiggs","cats4")) ;
	myWS.var("mHiggs_cat0")->setVal(mHiggs_cat0Val);
	myWS.var("mHiggs_cat1")->setVal(mHiggs_cat1Val);
	myWS.var("mHiggs_cat2")->setVal(mHiggs_cat2Val);
	myWS.var("mHiggs_cat3")->setVal(mHiggs_cat3Val);

	myWS.var("wHiggs_cat0")->setVal(wHiggs_cat0Val);
	myWS.var("wHiggs_cat1")->setVal(wHiggs_cat1Val);
	myWS.var("wHiggs_cat2")->setVal(wHiggs_cat2Val);
	myWS.var("wHiggs_cat3")->setVal(wHiggs_cat3Val);

	myWS.var("nbGauss_cat0")->setVal(nGauss_cat0);
	myWS.var("nbGauss_cat1")->setVal(nGauss_cat1);
	myWS.var("nbGauss_cat2")->setVal(nGauss_cat2);
	myWS.var("nbGauss_cat3")->setVal(nGauss_cat3);

	myWS.var("mHiggs2_cat0")->setVal(mHiggs2_cat0Val);
	myWS.var("mHiggs2_cat1")->setVal(mHiggs2_cat1Val);
	myWS.var("mHiggs2_cat2")->setVal(mHiggs2_cat2Val);
	myWS.var("mHiggs2_cat3")->setVal(mHiggs2_cat3Val);
	
	myWS.var("wHiggs2_cat0")->setVal(wHiggs2_cat0Val);
	myWS.var("wHiggs2_cat1")->setVal(wHiggs2_cat1Val);
	myWS.var("wHiggs2_cat2")->setVal(wHiggs2_cat2Val);
	myWS.var("wHiggs2_cat3")->setVal(wHiggs2_cat3Val);
	
	myWS.var("nbGauss2_cat0")->setVal(nGauss2_cat0);
	myWS.var("nbGauss2_cat1")->setVal(nGauss2_cat1);
	myWS.var("nbGauss2_cat2")->setVal(nGauss2_cat2);
	myWS.var("nbGauss2_cat3")->setVal(nGauss2_cat3);
	
	myWS.var("tau_cat0")->setVal(tau_cat0Val);
	myWS.var("tau_cat1")->setVal(tau_cat1Val);
	myWS.var("tau_cat2")->setVal(tau_cat2Val);
	myWS.var("tau_cat3")->setVal(tau_cat3Val);

	myWS.var("nbExp_cat0")->setVal(nbkg_cat0Val);
	myWS.var("nbExp_cat1")->setVal(nbkg_cat1Val);
	myWS.var("nbExp_cat2")->setVal(nbkg_cat2Val);
	myWS.var("nbExp_cat3")->setVal(nbkg_cat3Val);

	myWS.var("tau2_cat0")->setVal(tau2_cat0Val);
	myWS.var("tau2_cat1")->setVal(tau2_cat1Val);
	myWS.var("tau2_cat2")->setVal(tau2_cat2Val);
	myWS.var("tau2_cat3")->setVal(tau2_cat3Val);
	
	myWS.var("nbExp2_cat0")->setVal(nbkg2_cat0Val);
	myWS.var("nbExp2_cat1")->setVal(nbkg2_cat1Val);
	myWS.var("nbExp2_cat2")->setVal(nbkg2_cat2Val);
	myWS.var("nbExp2_cat3")->setVal(nbkg2_cat3Val);
	
	myWS.var("bkg_yield_cat0")->setVal(nbBg_cat0);
	myWS.var("bkg_yield_cat1")->setVal(nbBg_cat1);
	myWS.var("bkg_yield_cat2")->setVal(nbBg_cat2);
	myWS.var("bkg_yield_cat3")->setVal(nbBg_cat3);

	myWS.var("theNbHiggs_cat0")->setVal(nbHiggs_cat0);
	myWS.var("theNbHiggs_cat1")->setVal(nbHiggs_cat1);
	myWS.var("theNbHiggs_cat2")->setVal(nbHiggs_cat2);
	myWS.var("theNbHiggs_cat3")->setVal(nbHiggs_cat3);
	
	/*myWS.import(RooArgSet(bkg_yield_prior,cats4)) ;
	RooSimultaneous* bkg_yield_prior_sim = sct.build("tbkg_yield_prior_sim","bkg_yield_prior",SplitParam("bkg_yield_syst,bkg_yield_width","cats4")) ;
	for (int i = 0 ; i < nbOfCats ; i++){
		myWS.var(Form("bkg_yield_syst_cat%i",i))->setVal(nbBg[i]);
		myWS.var(Form("bkg_yield_width_cat%i",i))->setVal(nbBg[i]*5/100);
	}*/
	
/*myWS.var("bkg_yield_syst_cat0")->setVal(nbBg_cat0);		
	myWS.var("bkg_yield_syst_cat1")->setVal(nbBg_cat1);		
	myWS.var("bkg_yield_syst_cat2")->setVal(nbBg_cat2);		
	myWS.var("bkg_yield_syst_cat3")->setVal(nbBg_cat3);		
	
	myWS.var("bkg_yield_width_cat0")->setVal(nbBg_cat0*5/100);		
	myWS.var("bkg_yield_width_cat1")->setVal(nbBg_cat1*5/100);		
	myWS.var("bkg_yield_width_cat2")->setVal(nbBg_cat2*5/100);		
	myWS.var("bkg_yield_width_cat3")->setVal(nbBg_cat3*5/100);	*/
	
	RooArgSet nuisance_parameters(*myWS.var("bkg_yield_cat0"), *myWS.var("bkg_yield_cat1"),*myWS.var("bkg_yield_cat2"),*myWS.var("bkg_yield_cat3"));
	RooGaussian bkg_yield_prior_cat0("bkg_yield_prior_cat0","bkg_yield_prior_cat0",*myWS.var("bkg_yield_cat0"),RooConst(nbBg_cat0),RooConst(nbBg_cat0*5/100));
    RooGaussian bkg_yield_prior_cat1("bkg_yield_prior_cat1","bkg_yield_prior_cat1",*myWS.var("bkg_yield_cat1"),RooConst(nbBg_cat1),RooConst(nbBg_cat1*5/100));
    RooGaussian bkg_yield_prior_cat2("bkg_yield_prior_cat2","bkg_yield_prior_cat2",*myWS.var("bkg_yield_cat2"),RooConst(nbBg_cat2),RooConst(nbBg_cat2*5/100));
	RooGaussian bkg_yield_prior_cat3("bkg_yield_prior_cat3","bkg_yield_prior_cat3",*myWS.var("bkg_yield_cat3"),RooConst(nbBg_cat3),RooConst(nbBg_cat3*5/100));

	//RooProdPdf theTotSystPdf("theTotSystPdf","theTotSystPdf",RooArgSet(*myWS.pdf("bkg_yield_prior_cat0"),*myWS.pdf("bkg_yield_prior_cat1"), *myWS.pdf("bkg_yield_prior_cat2"),*myWS.pdf("bkg_yield_prior_cat3")));
	RooProdPdf theTotSystPdf("theTotSystPdf","theTotSystPdf",RooArgSet(bkg_yield_prior_cat0,bkg_yield_prior_cat1, bkg_yield_prior_cat2,bkg_yield_prior_cat3));

	myWS.import(theTotSystPdf);
	
	
	myWS.defineSet("obs",observables);
	myWS.defineSet("poi","r");
	myWS.defineSet("nuisances",nuisance_parameters);//, RenameConflictNodes("nui"));
	myWS.defineSet("cats4",cats4);
	
	//RooDataSet *data = tot_pdf_sim->generate(observables,RooFit::Extended(),Name("data"));  ExpectedData()
	myWS.var("r")->setVal(1.0);
	RooDataSet *data_cat0 = myWS.pdf("tot_pdf_cat0")->generate(observables,RooFit::Extended(),Name("datacat0"));
	RooDataSet *data_cat1 = myWS.pdf("tot_pdf_cat1")->generate(observables,RooFit::Extended(),Name("datacat1"));
	RooDataSet *data_cat2 = myWS.pdf("tot_pdf_cat2")->generate(observables,RooFit::Extended(),Name("datacat2"));
	RooDataSet *data_cat3 = myWS.pdf("tot_pdf_cat3")->generate(observables,RooFit::Extended(),Name("datacat3"));
	/*RooDataHist *data_cat0 = myWS.pdf("tot_pdf_cat0")->generateBinned(observables,RooFit::ExpectedData(),RooFit::Extended(),Name("datacat0"));
	RooDataHist *data_cat1 = myWS.pdf("tot_pdf_cat1")->generateBinned(observables,RooFit::ExpectedData(),RooFit::Extended(),Name("datacat1"));	
	RooDataHist *data_cat2 = myWS.pdf("tot_pdf_cat2")->generateBinned(observables,RooFit::ExpectedData(),RooFit::Extended(),Name("datacat2"));
	RooDataHist *data_cat3 = myWS.pdf("tot_pdf_cat3")->generateBinned(observables,RooFit::ExpectedData(),RooFit::Extended(),Name("datacat3"));*/
	

	RooDataSet *data_cat0 = new RooDataSet("datacat0","datacat0",chain0, mgg);
	RooDataSet *data_cat1 = new RooDataSet("datacat1","datacat1",chain1, mgg);
	RooDataSet *data_cat2 = new RooDataSet("datacat2","datacat2",chain2, mgg);
	RooDataSet *data_cat3 = new RooDataSet("datacat3","datacat3",chain3, mgg);


	/*RooAbsData *data_cat0 = (RooAbsData*) data.reduce("theCat==0") ;
	RooAbsData *data_cat1 = (RooAbsData*) data.reduce("theCat==1") ;
	RooAbsData *data_cat2 = (RooAbsData*) data.reduce("theCat==2") ;
	RooAbsData *data_cat3 = (RooAbsData*) data.reduce("theCat==3") ;*/



	myWS.import(*data_cat0);
	myWS.import(*data_cat1);
	myWS.import(*data_cat2);
	myWS.import(*data_cat3);
	
	//RooDataSet combData("combData","combined data",mgg,Index(cats4),Import("cat0",*data_cat0),Import("cat1",*data_cat1),Import("cat2",*data_cat2),Import("cat3",*data_cat3));
	//myWS.import(combData);
	myWS.Print("v") ;
	myWS.writeToFile("workSpace_M115_"+WSsuffixe+".root");
	myWS.var("r")->setVal(7.0);
	
	TCanvas *c0 = new TCanvas("c0","coucou",1200,1200);
	c0->Divide(2,2);
	c0->cd(1);
	RooPlot* frame1 = mgg.frame(Bins(30),Title("Physics sample")) ;
	data_cat0->plotOn(frame1);
	myWS.pdf("tot_pdf_cat0")->plotOn(frame1) ;
	frame1->Draw() ;

	c0->cd(2);
	RooPlot* frame2 = mgg.frame(Bins(30),Title("Physics sample")) ;
	data_cat1->plotOn(frame2);
	myWS.pdf("tot_pdf_cat1")->plotOn(frame2) ;
	frame2->Draw() ;
	
	c0->cd(3);
	RooPlot* frame3 = mgg.frame(Bins(30),Title("Physics sample")) ;
	data_cat2->plotOn(frame3);
	myWS.pdf("tot_pdf_cat2")->plotOn(frame3) ;
	frame3->Draw() ;
	
	c0->cd(4);
	RooPlot* frame4 = mgg.frame(Bins(30),Title("Physics sample")) ;
	data_cat3->plotOn(frame4);
	myWS.pdf("tot_pdf_cat3")->plotOn(frame4) ;
	frame4->Draw() ;
	
	
	c0->Print("gif/theModel.gif");

	

}

	
