#include "TString.h"
#include "TH1F.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>


TChain *chain = new TChain("diPhotons");

TH1F *th1fMaker(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TString nomEnPlus, TString cutEnPlus ){
	TString theCut;
/*	if (part==1) {theCut = "pho_isEB==1"; nom+="_EB";}
	else if (part==2) {theCut = "pho_isEE==1"; nom+="_EE";} 
	else {theCut = "1"; nom+="_All";}*/
	theCut += Cuts+"&&"+cutEnPlus;
	nom+="_"+nomEnPlus;
	TH1F *over = new TH1F("over","",(10000-Xmax), Xmax, 10000);
	TString plotOver = variable+">>over";
	chain->Draw(plotOver, theCut);
	int overFlow = over->Integral();
	delete over;
	float Lbin = (int) (Xmax - Xmin)/bining;
	float theLast = Xmax-Lbin;
	TH1F *h = new TH1F(nom,"",bining, Xmin, Xmax);
	TString toDraw = variable + ">>" + nom;
	chain->Draw(toDraw,theCut);
	h->Fill(theLast,overFlow);
	return h;
}

int doHistoGeneTroisCut(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TFile *theFile){
	cout << "processing " << nom << endl;
	TString localNom = "TwoReal";
	TString theLocalCutEnPlus = "((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))";
	TH1F *isrfsr = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);

        TString localNom = "OnePromptOneFake";
        TString theLocalCutEnPlus = "((((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&(!((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))||((!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	TH1F *prompt = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);

        TString localNom = "TwoFake";
        TString theLocalCutEnPlus = "(!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&(!(((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	TH1F *fake = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);

 	theFile->Write();

	delete isrfsr; delete prompt; delete fake;
}

int do_hist_2photons_MC(){
  TFile *myFile = new TFile("dipho_file.root","RECREATE");
   chain->Add("diphoton_QCD20.root");

TString ptHatCut = "event_eventPtHat>0&&";  
TString baseCut = "event_nVertex==1&&(abs(pholead_eta)<=2.5)&&(!((abs(pholead_eta)>1.4442)&&(abs(pholead_eta)<1.566)))&&(abs(pholead_eta)<=2.5)&&(!((abs(pholead_eta)>1.4442)&&(abs(pholead_eta)<1.566)))&&pholead_hoe<0.05&&photrail_hoe<0.05&&pholead_pt>30&&photrail_pt>30";
TString isoLead = "(pholead_TrackerIso<2&&pholead_EcalIso<4.2&&pholead_HcalIso<2.2)";
TString isoTrail = "(photrail_TrackerIso<2&&photrail_EcalIso<4.2&&photrail_HcalIso<2.2)";
TString cutIeta = "(((pholead_isEB==1)&&(pholead_sigieta<0.01))||((pholead_isEE==1)&&(pholead_sigieta<0.03)))&&(((photrail_isEB==1)&&(photrail_sigieta<0.01))||((photrail_isEE==1)&&(photrail_sigieta<0.03)))";
TString pixelCut = "pholead_HasPixSeed==0&&photrail_HasPixSeed==0";
TString combisoLead = "((pholead_TrackerIso+pholead_EcalIso+pholead_HcalIso)<3.0)";
TString combisoTrail = "((photrail_TrackerIso+photrail_EcalIso+photrail_HcalIso)<3.0)";

baseCut = ptHatCut + baseCut;

TString noIetaPixel = cutIeta + "&&" + pixelCut;

doHistoGeneTroisCut("leadEt_noIso","pholead_pt",0,150,30,noIetaPixel, myFile);
doHistoGeneTroisCut("trailEt_noIso","photrail_pt",0,150,30,noIetaPixel, myFile);

doHistoGeneTroisCut("isoTrackerLead","pholead_TrackerIso",-1,20,42,noIetaPixel, myFile);
doHistoGeneTroisCut("isoTrackerTrail","photrail_TrackerIso",-1,20,42,noIetaPixel, myFile);

doHistoGeneTroisCut("isoEcalLead","pholead_EcalIso",-2,16,36,noIetaPixel, myFile);
doHistoGeneTroisCut("isoEcalTrail","photrail_EcalIso",-2,16,36,noIetaPixel, myFile);

doHistoGeneTroisCut("isoHcalLead","pholead_HcalIso",-1,10,22,noIetaPixel, myFile);
doHistoGeneTroisCut("isoHcalTrail","photrail_HcalIso",-1,10,22,noIetaPixel, myFile);

doHistoGeneTroisCut("sumIsoLead","(pholead_TrackerIso+pholead_EcalIso+pholead_HcalIso)",-3,40,43,noIetaPixel, myFile);
doHistoGeneTroisCut("sumIsoTrail","(photrail_TrackerIso+photrail_EcalIso+photrail_HcalIso)",-3,40,43,noIetaPixel, myFile);

doHistoGeneTroisCut("cosThetaStarLead_noIso","dipho_costhetastar",-1,1,20,noIetaPixel, myFile);

doHistoGeneTroisCut("etaStar_noIso","dipho_etastar",-3,3,24,noIetaPixel, myFile);

doHistoGeneTroisCut("qt_noIso","dipho_qt",0,150,30,noIetaPixel, myFile);

doHistoGeneTroisCut("ql_noIso","dipho_ql",-750,750,30,noIetaPixel, myFile);


////////////////////////////////////////////////////////////////////// after isolation cuts
TString ietaIso = cutIeta +"&&" + isoLead+ "&&" + isoTrail + "&&" + pixelCut;
doHistoGeneTroisCut("leadEt_iso","pholead_pt",0,150,30,ietaIso, myFile);
doHistoGeneTroisCut("trailEt_iso","photrail_pt",0,150,30,ietaIso, myFile);
doHistoGeneTroisCut("cosThetaStarLead_iso","dipho_costhetastar",-1,1,20,ietaIso, myFile);

doHistoGeneTroisCut("etaStar_iso","dipho_etastar",-3,3,24,ietaIso, myFile);
doHistoGeneTroisCut("qt_iso","dipho_qt",0,150,30,ietaIso, myFile);
doHistoGeneTroisCut("ql_iso","dipho_ql",-750,750,30,ietaIso, myFile);


////////////////////////////////////////////////////////////////////// after isolation cuts
TString ietaCombIso = cutIeta +"&&" + isoLead+ "&&" + isoTrail + "&&" + pixelCut;
doHistoGeneTroisCut("leadEt_combIso","pholead_pt",0,150,30,ietaCombIso, myFile);
doHistoGeneTroisCut("trailEt_combIso","photrail_pt",0,150,30,ietaCombIso, myFile);
doHistoGeneTroisCut("cosThetaStarLead_combIso","dipho_costhetastar",-1,1,20,ietaCombIso, myFile);

doHistoGeneTroisCut("etaStar_combIso","dipho_etastar",-3,3,24,ietaCombIso, myFile);
doHistoGeneTroisCut("qt_combIso","dipho_qt",0,150,30,ietaCombIso, myFile);
doHistoGeneTroisCut("ql_combIso","dipho_ql",-750,750,30,ietaCombIso, myFile);



}
