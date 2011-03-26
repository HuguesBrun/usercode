#include "TString.h"
#include "TH1F.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>
TString name = "data";

TChain *chain = new TChain("photons");

int printCutFlow(TString *cuts, int nbCuts, TString nameTag, TString localCut){
	
	TH1F *h = new TH1F("h","",100,-3,3);
	chain->Draw("pho_eta>>h");
	cout << "<" << nameTag << ">";
	cout << "<cut0>" << h->Integral() << "</cut0>";
	TString theCuts = "";
	theCuts = localCut+"&&";
	for (int i = 0 ; i < nbCuts; i++){
		theCuts+=cuts[i];
		chain->Draw("pho_eta>>h",theCuts);
		cout << "<cut" << i+1 << ">" << h->Integral() << "</cut" << i+1 << ">";
		theCuts+="&&";
	}
	delete h;
	cout << "</" << nameTag << ">" << endl;
}

TH1F *th1fMaker(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, int part, TString nomEnPlus, TString cutEnPlus ){
	TString theCut;
	if (part==1) {theCut = "pho_isEB==1"; nom+="_EB";}
	else if (part==2) {theCut = "pho_isEE==1"; nom+="_EE";} 
	else {theCut = "1"; nom+="_All";}
	theCut += "&&"+Cuts+"&&"+cutEnPlus;
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

int doHistoGeneTroisCut(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, int part, TFile *theFile){
        cout << "processing the plot " << nom << endl;
	TH1F *data = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,part, "1", "1");
 	theFile->Write();
	delete data;
}

int doHistoGeneTroisCut_SC(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, int part, TFile *theFile){
        cout << "processing the plot " << nom << endl;
	TH1F *data = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,part,"1", "1");
 	theFile->Write();
	delete data;
}

TH1F *doNminus1(TString nom, TString variable, float Xmin, float Xmax, float bining, TString* cuts, int nbCut, int cutNb,int part, TString nomEnPlus, TString cutEnPlus){
	TString theNminus1Cut;
	if (part==1) {theNminus1Cut = "pho_isEB==1"; nom+="_EB";}
	else if (part==2) {theNminus1Cut = "pho_isEE==1"; nom+="_EE";} 
	else {theNminus1Cut = "1"; nom+="_All";}
	theNminus1Cut += "&&"+cutEnPlus;
	nom+="_"+nomEnPlus;
	TH1F *over = new TH1F("over","",(1000-Xmax), Xmax, 1000);
	TH1F *h = new TH1F(nom,"",bining, Xmin, Xmax);
	for (int i = 0 ; i < nbCut ; i++){
		if (i == cutNb) continue;
		theNminus1Cut += "&&";
		theNminus1Cut += cuts[i];
	}
	TString plotOver = variable+">>over";
	chain->Draw(plotOver, theNminus1Cut);
	int overFlow = over->Integral();
	delete over;
	float Lbin = (int) (Xmax - Xmin)/bining;
	float theLast = Xmax-Lbin;
	TString toDraw = variable + ">>" + nom;
	chain->Draw(toDraw,theNminus1Cut);
	h->Fill(theLast,overFlow);
	return h;
}

int doHistoGeneTroisCutNminus1(TString nom, TString variable, float Xmin, float Xmax, float bining, TString* cuts, int nbCut, int cutNb,int part, TFile *theFile){
        cout << "processing the plot " << nom << endl;
	TH1F *data = (TH1F*) doNminus1(nom,variable,Xmin,Xmax,bining, cuts, nbCut, cutNb, part,"1","1");
 	theFile->Write();

	delete data;
}


int do_hist(){
//int main(){
TFile *myFile = new TFile("histo_file.root","RECREATE");

chain->Add("/sps/cms/hbrun/miniTree38X/part3/output_1.root");

int cutOffSet = 1;
int nbCuts = 9;
TString cuts[9] = {"pho_HLT_bit1==1","isAspike==0&&(abs(pho_eta)<=2.5)&&(!((abs(pho_eta)>1.4442)&&(abs(pho_eta)<1.566)))","pho_hoe<0.05","pho_et>30","pho_hasPixelSeed==0","pho_IsoHollowTrkCone<2","pho_IsoEcalRechit<4.2","pho_IsoHcalRechit<2.2","((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))"};
////////////////////   Modif pour les cut en pt hat

int nbCutsIso =7;
TString cutsIso[7] = {"pho_HLT_bit1==1","isAspike==0&&(abs(pho_eta)<=2.5)&&(!((abs(pho_eta)>1.4442)&&(abs(pho_eta)<1.566)))","pho_hasPixelSeed==0","pho_et>30","pho_IsoHollowTrkCone<2&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2","pho_hoe<0.05","((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))"};

		printCutFlow(cuts, nbCuts, name, "1");


TString AlltheCuts = "";
for (int i = 0 ; i < nbCuts; i++){
         AlltheCuts+=cuts[i];
         AlltheCuts+="&&";
 }
AlltheCuts+="1";

	doHistoGeneTroisCutNminus1("pixelSeed","pho_hasPixelSeed",0,2,2,cuts, nbCuts, 3+cutOffSet,1, myFile );
	doHistoGeneTroisCutNminus1("pixelSeed","pho_hasPixelSeed",0,2,2,cuts, nbCuts, 3+cutOffSet,2, myFile );
	
	doHistoGeneTroisCutNminus1("sumIso","(pho_IsoHollowTrkCone+pho_IsoEcalRechit+pho_IsoHcalRechit)",-3,40,43,cutsIso, nbCutsIso, 3+cutOffSet,1, myFile);
	doHistoGeneTroisCutNminus1("sumIso","(pho_IsoHollowTrkCone+pho_IsoEcalRechit+pho_IsoHcalRechit)",-3,40,43,cutsIso, nbCutsIso, 3+cutOffSet,2, myFile);
	
	doHistoGeneTroisCutNminus1("isoTrack","pho_IsoHollowTrkCone",-1,20,42,cuts, nbCuts, 4+cutOffSet,1, myFile);
	doHistoGeneTroisCutNminus1("isoTrack","pho_IsoHollowTrkCone",-1,20,42,cuts, nbCuts, 4+cutOffSet,2, myFile);

	doHistoGeneTroisCutNminus1("isoEcal","pho_IsoEcalRechit",-2,16,36,cuts, nbCuts, 5+cutOffSet,1, myFile);
	doHistoGeneTroisCutNminus1("isoEcal","pho_IsoEcalRechit",-2,16,36,cuts, nbCuts, 5+cutOffSet,2, myFile);
	
	doHistoGeneTroisCutNminus1("isoHcal","pho_IsoHcalRechit",-1,10,22,cuts, nbCuts, 6+cutOffSet,1, myFile);
	doHistoGeneTroisCutNminus1("isoHcal","pho_IsoHcalRechit",-1,10,22,cuts, nbCuts, 6+cutOffSet,2, myFile);
	
	doHistoGeneTroisCutNminus1("sigmaIeta","pho_sigmaIetaIeta",0,0.04,50,cuts, nbCuts, 7+cutOffSet,1,  myFile);
	doHistoGeneTroisCutNminus1("sigmaIeta","pho_sigmaIetaIeta",0,0.01,50,cuts, nbCuts, 7+cutOffSet,2,  myFile);
	
	
	doHistoGeneTroisCut("Et","pho_et",0,250,50,AlltheCuts,1, myFile);
	doHistoGeneTroisCut("Et","pho_et",0,250,50,AlltheCuts,2, myFile);
	
	doHistoGeneTroisCut("photon_r9","pho_r9",0,1.2,40,AlltheCuts,1, myFile);
	doHistoGeneTroisCut("photon_r9","pho_r9",0,1.2,40,AlltheCuts,2, myFile);	
	
	doHistoGeneTroisCut("photon_eta","pho_eta",-3.1,3.1,31,AlltheCuts,0,myFile);
	
	
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int nbCutsSC = 4;
TString cutsSC[4] = {"pho_HLT_bit1==1","(abs(pho_SCeta)<=2.5)&&(!((abs(pho_SCeta)>1.4442)&&(abs(pho_SCeta)<1.566)))","pho_hoe<0.05","pho_SCEt>20"};

printCutFlow(cutsSC, nbCutsSC);

TString AlltheCutsSC = "";
for (int i = 0 ; i < nbCutsSC; i++){
         AlltheCutsSC+=cutsSC[i];
         AlltheCutsSC+="&&";
 }
AlltheCutsSC+="1";

	    doHistoGeneTroisCut_SC("SC_eta","pho_SCeta",-3.1,3.1,32,AlltheCutsSC,0, myFile);

	    doHistoGeneTroisCut_SC("SC_Et","pho_SCEt",0,250,50,AlltheCutsSC,1, myFile);
	    doHistoGeneTroisCut_SC("SC_Et","pho_SCEt",0,250,50,AlltheCutsSC,2, myFile);
	    
	    doHistoGeneTroisCut_SC("SC_r9","pho_SCr9",0,1.3,52,AlltheCutsSC,1, myFile);
	    doHistoGeneTroisCut_SC("SC_r9","pho_SCr9",0,1.3,52,AlltheCutsSC,2, myFile);

	    doHistoGeneTroisCut_SC("SC_nXtal","pho_SCnXtal",0,500,500,AlltheCutsSC,1, myFile);
	    doHistoGeneTroisCut_SC("SC_nXtal","pho_SCnXtal",0,500,500,AlltheCutsSC,2, myFile);
	    
	    doHistoGeneTroisCut_SC("SC_BC","pho_SCnbBC",0,15,15,AlltheCutsSC,1, myFile);
	    doHistoGeneTroisCut_SC("SC_BC","pho_SCnbBC",0,15,15,AlltheCutsSC,2, myFile);
	    
	    doHistoGeneTroisCut_SC("SC_br","pho_SCbr",0,15,75,AlltheCutsSC,1, myFile);
	    doHistoGeneTroisCut_SC("SC_br","pho_SCbr",0,15,75,AlltheCutsSC,2, myFile);
}
