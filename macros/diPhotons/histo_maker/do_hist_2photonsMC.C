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
TString baseCut = "(abs(pholead_eta)<=2.5)&&(!((abs(pholead_eta)>1.4442)&&(abs(pholead_eta)<1.566)))&&(abs(photrail_eta)<=2.5)&&(!((abs(photrail_eta)>1.4442)&&(abs(photrail_eta)<1.566)))&&pholead_hoe<0.05&&photrail_hoe<0.05&&pholead_pt>30&&photrail_pt>30";
TString isoLead = "(pholead_TrackerIso<2&&pholead_EcalIso<4.2&&pholead_HcalIso<2.2)";
TString isoTrail = "(photrail_TrackerIso<2&&photrail_EcalIso<4.2&&photrail_HcalIso<2.2)";
TString cutIeta = "(((pholead_isEB==1)&&(pholead_sigieta<0.01))||((pholead_isEE==1)&&(pholead_sigieta<0.03)))&&(((photrail_isEB==1)&&(photrail_sigieta<0.01))||((photrail_isEE==1)&&(photrail_sigieta<0.03)))";
TString pixelCut = "pholead_HasPixSeed==0&&photrail_HasPixSeed==0";
TString combisoLead = "((pholead_TrackerIso+pholead_EcalIso+pholead_HcalIso)<3.0)";
TString combisoTrail = "((photrail_TrackerIso+photrail_EcalIso+photrail_HcalIso)<3.0)";

baseCut = ptHatCut + baseCut;

TString noIetaPixel = baseCut + "&&" + pixelCut;
doHistoGeneTroisCut("noIetaPixel","dipho_mgg",0,510,34,noIetaPixel, myFile);

TString baseCutPixel = baseCut+"&&"+pixelCut+"&&" +cutIeta;
doHistoGeneTroisCut("baseCutPixel","dipho_mgg",0,510,34,baseCutPixel, myFile);

TString isoLeadPixel = baseCut+"&&"+pixelCut+"&&"+isoLead+"&&" +cutIeta;
doHistoGeneTroisCut("leadSelecPixel","dipho_mgg",0,519,34,isoLeadPixel, myFile);

TString isoTrailPixel = baseCut+"&&"+pixelCut+"&&"+isoTrail+"&&" +cutIeta;
doHistoGeneTroisCut("trailSelecPixel","dipho_mgg",0,510,34,isoTrailPixel, myFile);

TString isoOnePixel = baseCut+"&&"+pixelCut+"&&("+isoLead + "||" +isoTrail+")"+"&&" +cutIeta;
doHistoGeneTroisCut("oneSelecPixel","dipho_mgg",0,510,34,isoOnePixel, myFile);               

TString isoBothPixel = baseCut+"&&"+pixelCut+"&&("+isoLead + "&&" +isoTrail+")"+"&&" +cutIeta;
doHistoGeneTroisCut("bothSelecPixel","dipho_mgg",0,510,34,isoBothPixel, myFile);              

TString combisoLeadPixel  =baseCut + "&&" + pixelCut + "&&" + combisoLead+"&&" +cutIeta;
doHistoGeneTroisCut("leadCombIsoPixel","dipho_mgg",0,510,34,combisoLeadPixel, myFile);  

TString combisoTrailPixel  =baseCut + "&&" + pixelCut + "&&" + combisoTrail+"&&" +cutIeta;
doHistoGeneTroisCut("trailCombIsoPixel","dipho_mgg",0,510,34,combisoTrailPixel, myFile);  

TString combisoOnePixel  =baseCut + "&&" + pixelCut + "&&(" + combisoLead+"||" + combisoTrail + ")&&" +cutIeta;
doHistoGeneTroisCut("oneCombIsoPixel","dipho_mgg",0,510,34,combisoOnePixel, myFile);                           

TString combisoBothPixel  =baseCut + "&&" + pixelCut + "&&(" + combisoLead+"&&" + combisoTrail + ")&&" +cutIeta;
doHistoGeneTroisCut("bothCombIsoPixel","dipho_mgg",0,510,34,combisoBothPixel, myFile);                          

///////////////////////////////////////////////////Now with no pixel !!

TString noIeta = baseCut;
doHistoGeneTroisCut("noIetaNoPixel","dipho_mgg",0,510,102,baseCut, myFile);

TString baseCut = baseCut+"&&" +cutIeta;
doHistoGeneTroisCut("baseCutNoPixel","dipho_mgg",0,510,102,baseCut, myFile);

TString isoLead = baseCut+"&&"+isoLead+"&&" +cutIeta;
doHistoGeneTroisCut("leadSelecNoPixel","dipho_mgg",0,510,102,isoLead, myFile);

TString isoTrail = baseCut+"&&"+isoTrail+"&&" +cutIeta;
doHistoGeneTroisCut("trailSelecNoPixel","dipho_mgg",0,510,102,isoTrail, myFile);

TString isoOne = baseCut+"&&("+isoLead + "||" +isoTrail+")"+"&&" +cutIeta;
doHistoGeneTroisCut("oneSelecNoPixel","dipho_mgg",0,510,102,isoOne, myFile);

TString isoBoth = baseCut+"&&("+isoLead + "&&" +isoTrail+")"+"&&" +cutIeta;
doHistoGeneTroisCut("bothSelecNoPixel","dipho_mgg",0,510,102,isoBoth, myFile);

TString combisoLead =baseCut + "&&" + combisoLead+"&&" +cutIeta;
doHistoGeneTroisCut("leadCombIsoNoPixel","dipho_mgg",0,510,102,combisoLead, myFile);

TString combisoTrail  =baseCut + "&&" +  combisoTrail+"&&" +cutIeta;
doHistoGeneTroisCut("trailCombIsoNoPixel","dipho_mgg",0,510,102,combisoTrail, myFile);

TString combisoOne  =baseCut + "&&(" + combisoLead+"||" + combisoTrail + ")&&" +cutIeta;
doHistoGeneTroisCut("oneCombIsoNoPixel","dipho_mgg",0,510,102,combisoOne, myFile);

TString combisoBoth  =baseCut + "&&(" + combisoLead+"&&" + combisoTrail + ")&&" +cutIeta;
doHistoGeneTroisCut("bothCombIsoNoPixel","dipho_mgg",0,510,102,combisoBoth, myFile);

TString combisoBothWithPixel  =baseCut + "&&(" + combisoLead+"&&" + combisoTrail + ")&&" +cutIeta + "&&pholead_HasPixSeed==1&&photrail_HasPixSeed==1";
doHistoGeneTroisCut("bothCombIsoWithPixel","dipho_mgg",0,510,102,combisoBothWithPixel, myFile);


// do a plot for the Z fit !!!

doHistoGeneTroisCut("bothCombIsoWithPixelZoom","dipho_mgg",70,110,30,combisoBothWithPixel, myFile);



}
