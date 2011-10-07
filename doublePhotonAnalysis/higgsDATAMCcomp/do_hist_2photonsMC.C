#include "TString.h"
#include "TH1F.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>

TString name = "mc";

TString FloatToString(float number){
	ostringstream oss;
	oss << number;
	return oss.str();
}


TString IntToString(int number){
	ostringstream oss;
	oss << number;
	return oss.str();
}



TChain *chain = new TChain("diPhotons");

float giveTheWeight(int nbVertex){
	if (nbVertex > 25) return 0;
	float theWeight[25]  = {0, 0.186814, 1.06775, 1.73441, 1.98779, 2.20497, 1.92209, 1.68128, 1.22481, 0.884089, 0.64283, 0.47647, 0.276678, 0.169805, 0.115875, 0.0630132, 0.0376695, 0.0225988, 0.0152533, 0.0094274, 0.00500494, 0.00290434, 0.00127482, 0.000728612, 0.00077654};
        return theWeight[nbVertex];
}


int printCutFlow(TString *cuts, int nbCuts, TString nameTag, TString localCut){
	
	TH1F *h = new TH1F("h","",100,-3,3);
	chain->Draw("pholead_eta>>h","giveTheWeight(event_nVertex)");
	cout << "<" << nameTag << ">";
	cout << "<cut0>" << h->Integral() << "</cut0>";
	TString theCuts = "";
	theCuts = localCut+"&&";
	for (int i = 0 ; i < nbCuts; i++){
		theCuts+=cuts[i];
		TString theCutsWeighted = "(" + theCuts + ")";//*giveTheWeight(event_nVertex)";
		chain->Draw("pholead_eta>>h",theCutsWeighted);
		cout << "<cut" << i+1 << ">" << h->Integral() << "</cut" << i+1 << ">";
		theCuts+="&&";
	}
	delete h;
	cout << "</" << nameTag << ">";
}



int printAllCutFlow(TString *cuts, int nbCuts){
	cout << "<cutFlow_"+name+">" << endl;
	printCutFlow(cuts,nbCuts,"all","1");
	for (int i = 90; i < 150 ; i=i+5){
		int lowRange = i - 5;
		int highRange = i + 5;
		TString lowPart = "dipho_mgg>"+FloatToString(lowRange);
		TString highPart = "dipho_mgg<"+FloatToString(highRange);
		TString theCond = lowPart + "&&" + highPart;
		TString theName = "masse"+FloatToString(i);
		printCutFlow(cuts,nbCuts,theName,theCond);
	}
	cout << "</cutFlow_"+name+">" << endl;
}



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


TH1F *th1fMakerWeight(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TString nomEnPlus, TString cutEnPlus){

        /*if (part==1) {theCut = thePhoton+"_isEB==1"; nom+="_EB";}
        else if (part==2) {theCut = thePhoton+"_isEE==1"; nom+="_EE";}
        else {theCut = "1"; nom+="_All";}*/

        nom+="_"+nomEnPlus;
        TH1F *h = new TH1F(nom,"",bining, Xmin, Xmax);
        TString OriCuts = Cuts;
        for (int i = 0; i<25; i++){
           Cuts = Cuts + "&& event_nGenVertex == " + IntToString(i);
           cout << "nGenVertex == " << IntToString(i)<<endl;
           cout << Cuts <<endl;
           TH1F *histoBin = th1fMaker(nom, variable, Xmin, Xmax, bining, Cuts, nomEnPlus, cutEnPlus);
           h->Add(histoBin, giveTheWeight(i));
           delete histoBin;
           Cuts = OriCuts;

        }
        return h;
}


int doHistoGeneTroisCut(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TFile *theFile){
	cout << "processing " << nom << endl;
	TString localNom = "TwoReal";
	TString theLocalCutEnPlus = "((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))";
	//TH1F *isrfsr = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	TH1F *isrfsr = (TH1F*) th1fMakerWeight(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);

        TString localNom = "OnePromptOneFake";
        TString theLocalCutEnPlus = "((((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&(!((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))||((!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	//TH1F *prompt = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	TH1F *prompt = (TH1F*) th1fMakerWeight(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);

        TString localNom = "TwoFake";
        TString theLocalCutEnPlus = "(!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&(!(((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	//TH1F *fake = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	TH1F *fake = (TH1F*) th1fMakerWeight(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);

 	theFile->Write();

	delete isrfsr; delete prompt; delete fake;
}


TH1F *th1fMakerNoOverF(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TString nomEnPlus, TString cutEnPlus ){
	TString theCut;
	/*	if (part==1) {theCut = "pho_isEB==1"; nom+="_EB";}
	 else if (part==2) {theCut = "pho_isEE==1"; nom+="_EE";} 
	 else {theCut = "1"; nom+="_All";}*/
	theCut += Cuts+"&&"+cutEnPlus;
	nom+="_"+nomEnPlus;
	TString theCutsWeighted = "(" + theCut + ")";//*giveTheWeight(event_nVertex)";
	TH1F *h = new TH1F(nom,"",bining, Xmin, Xmax);
	TString toDraw = variable + ">>" + nom;
	chain->Draw(toDraw,theCutsWeighted);
	return h;
}


TH1F *th1fMakerNoOverFWeight(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TString nomEnPlus, TString cutEnPlus ){
	/*	if (part==1) {theCut = "pho_isEB==1"; nom+="_EB";}
	 else if (part==2) {theCut = "pho_isEE==1"; nom+="_EE";} 
	 else {theCut = "1"; nom+="_All";}*/
	nom+="_"+nomEnPlus;
 
        TH1F *h = new TH1F(nom,"",bining, Xmin, Xmax);
        TString OriCuts = Cuts;
        for (int i = 0; i<25; i++){
           Cuts = Cuts + "&& event_nGenVertex == " + IntToString(i);
           cout << "nGenVertex == " << IntToString(i)<<endl;
           cout << Cuts <<endl;
           TH1F *histoBin = th1fMakerNoOverF(nom, variable, Xmin, Xmax, bining, Cuts, nomEnPlus, cutEnPlus);
           h->Add(histoBin, giveTheWeight(i));
           delete histoBin;
           Cuts = OriCuts;
        }

        return h;
}


int doHistoGeneTroisCutNoOverF(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, TFile *theFile){
	cout << "processing " << nom << endl;
	TString localNom = "TwoReal";
	TString theLocalCutEnPlus = "((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))";
	//TH1F *isrfsr = (TH1F*) th1fMakerNoOverF(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	TH1F *isrfsr = (TH1F*) th1fMakerNoOverFWeight(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	
	TString localNom = "OnePromptOneFake";
	TString theLocalCutEnPlus = "((((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&(!((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))||((!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	//TH1F *prompt = (TH1F*) th1fMakerNoOverF(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	TH1F *prompt = (TH1F*) th1fMakerNoOverFWeight(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	
	TString localNom = "TwoFake";
	TString theLocalCutEnPlus = "(!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&(!(((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	//TH1F *fake = (TH1F*) th1fMakerNoOverF(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	TH1F *fake = (TH1F*) th1fMakerNoOverFWeight(nom,variable,Xmin,Xmax,bining,Cuts,localNom, theLocalCutEnPlus);
	
 	theFile->Write();
	
	delete isrfsr; delete prompt; delete fake;
}



TH1F *doNminus1(TString nom, TString variable, float Xmin, float Xmax, float bining, TString* cuts, int nbCut, int cutNb, TString nomEnPlus, TString cutEnPlus){
	TString theNminus1Cut= "1";
/*	if (part==1) {theNminus1Cut = "pho_isEB==1"; nom+="_EB";}
	else if (part==2) {theNminus1Cut = "pho_isEE==1"; nom+="_EE";} 
	else {theNminus1Cut = "1"; nom+="_All";}*/
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
	TString theNminus1CutWeighted = "(" + theNminus1Cut + ")";//*giveTheWeight(event_nVertex)";
	chain->Draw(toDraw,theNminus1CutWeighted);
	h->Fill(theLast,overFlow);
	return h;
}

int doHistoGeneTroisCutNminus1(TString nom, TString variable, float Xmin, float Xmax, float bining, TString* cuts, int nbCut, int cutNb, TFile *theFile){
        cout << "processing the plot " << nom << endl;
    	TString localNom = "TwoReal";
	TString theLocalCutEnPlus = "((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))";
	TH1F *isrfsr = (TH1F*) doNminus1(nom,variable,Xmin,Xmax,bining,cuts,nbCut,cutNb,localNom, theLocalCutEnPlus);

        TString localNom = "OnePromptOneFake";
        TString theLocalCutEnPlus = "((((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&(!((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))||((!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	TH1F *prompt = (TH1F*) doNminus1(nom,variable,Xmin,Xmax,bining,cuts,nbCut,cutNb,localNom, theLocalCutEnPlus);

        TString localNom = "TwoFake";
        TString theLocalCutEnPlus = "(!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&(!(((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	TH1F *fake = (TH1F*) doNminus1(nom,variable,Xmin,Xmax,bining,cuts,nbCut,cutNb,localNom, theLocalCutEnPlus);

 	theFile->Write();

	delete isrfsr; delete prompt; delete fake;    
        
}


int do_hist_2photonsMC(){
  TFile *myFile = new TFile("dipho_file.root","RECREATE");
   chain->Add("diphoton_DiPhoBox25to250.root");

TString theCutAvant= "";


int theSelecType=1;
TString prefixe;
int nbCuts = 17;
TString selection[17];
if (theSelecType==0){
        prefixe = "loose1";
        int nbCuts = 17;
        selection[0] = "(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&dipho_mgg>85"+theCutAvant;
        selection[1] = 	"(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))";
        selection[2] = 	"(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))";
        selection[3] = 	"(((pholead_isEB==1)&&(pholead_sigieta<0.013))||((pholead_isEE==1)&&(pholead_sigieta<0.034)))";
        selection[4] = 	"(((photrail_isEB==1)&&(photrail_sigieta<0.013))||((photrail_isEE==1)&&(photrail_sigieta<0.034)))";
        selection[5] = 	"pholead_HasPixSeed==0";
        selection[6] = 	"photrail_HasPixSeed==0";
        selection[7] = 	"pholead_pt>40";
        selection[8] = 	"photrail_pt>30";
        selection[9] = 	"pholead_TrackerIso<2";
        selection[10] = 	"photrail_TrackerIso<2";
        selection[11] = 	"pholead_EcalIso<4.2";
        selection[12] = 	"photrail_EcalIso<4.2";
        selection[13] = 	"pholead_HcalIso<2.2";
        selection[14] = 	"photrail_HcalIso<2.2";
        selection[15] = 	"pholead_hoe<0.05";
        selection[16] = 	"photrail_hoe<0.05";
}
else if (theSelecType==1){
        prefixe = "CIC";
        int nbCuts = 7;
	selection[0] = "(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&dipho_mgg>85"+theCutAvant;
	selection[1] = "pholead_pt>40";
	selection[2] = "photrail_pt>30";
	selection[3] = "(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))";
 	selection[4] = "(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))";
        selection[5] = "pholead_isPassingCIC==1";
        selection[6] = "photrail_isPassingCIC==1";
}


	int nbVars = 29;
	TString varsToPlot[29]  = {
	        "dipho_mgg",
                "dipho_qt",
                "dipho_costhetastar",
                "dipho_etastar", 
                "pholead_eta",
                "photrail_eta",
                "pholead_pt",
                "photrail_pt",
                "pholead_r9",
                "photrail_r9",
                "pholead_NNshapeOutput",
                "photrail_NNshapeOutput",
                "dipho_minNNshape",
                "pholead_cEP",
                "photrail_cEP",
                "pholead_etawidth",
                "photrail_etawidth",
                "pholead_r19",
                "photrail_r19",
                "pholead_SCbr",
                "photrail_SCbr",
                "pholead_ratioSeed",
                "photrail_ratioSeed",
                "pholead_ratioS4",
                "photrail_ratioS4",
                "pholead_lambdaRatio",
                "photrail_lambdaRatio",
                "pholead_lamdbaDivCov",
                "photrail_lamdbaDivCov"
        };
	
	double bining[87] = {
		80,150,14, //Mgg
		0,150,30, // qT
		0,1,20, //cos(ThetaStar)
		-3,3,25, //etaStar
		-3.1,3.1,32, //photon eta
		-3.1,3.1,32,
		0,150,30, // photon pt
		0,150,30,
		0,1.2,40, // photon R9
		0,1.2,40,
		-0.2,1.2,28, // photon NN output
		-0.2,1.2,28,
		-0.2,1.2,28,
		-0.0005,0.0005,50, //cEP 
		-0.0005,0.0005,50,
		0.005,0.035,30,   //eta width
		0.005,0.035,30,
		0.2,0.9,35,  // r19
		0.2,0.9,35,
		0,8,40,     // brem
		0,8,40,
		0.1,0.9,40,  // Emax / Esc
		0.1,0.9,40,
		0.5,1,50,  // s4 ratio
		0.5,1,50,
		0,1,40,  // lamdba ratio
		0,1,40,
		1.2,2.1,45,  // lambda div cov 
		1.2,2.1,45
	};
	
	
	
	

	

//printAllCutFlow(selection, nbCuts);

TString allTheCuts = "1";
cout << "coucou = " << allTheCuts << endl;
for (int i = 0 ; i < nbCuts ; i++){
	allTheCuts += "&&" + selection[i];
	if (i<(nbCuts-1)) continue; 
	for (int j = 0 ; j < nbVars ; j++){
		TString nom = varsToPlot[j] + "_" + prefixe + Form("Cut%i",i);
		cout << allTheCuts << endl;
                //doHistoGeneTroisCut(Form("higgsTest_Cut%i",i),"dipho_mgg",90,150,60,allTheCuts,myFile);
		doHistoGeneTroisCut(nom,varsToPlot[j],bining[j*3],bining[j*3+1],bining[j*3+2],allTheCuts,myFile);
		cout << "apres" << endl;
	}
	
}
        doHistoGeneTroisCut("vertex","event_nVertex",0,25,25,allTheCuts,myFile);


       
	doHistoGeneTroisCutNoOverF("higgs1_"+prefixe,"dipho_mgg",70,150,80,allTheCuts,myFile);
	doHistoGeneTroisCutNoOverF("higgs2_"+prefixe,"dipho_mgg",80,150,35,allTheCuts,myFile);
	doHistoGeneTroisCutNoOverF("higgs3_"+prefixe,"dipho_mgg",90,150,20,allTheCuts,myFile);
	doHistoGeneTroisCutNoOverF("higgs4_"+prefixe,"dipho_mgg",90,150,15,allTheCuts,myFile);
	doHistoGeneTroisCutNoOverF("higgs5_"+prefixe,"dipho_mgg",80,150,14,allTheCuts,myFile);
        






/*	TString allTheCutsNN = allTheCuts + "&&dipho_minNNshape>0.7";
	doHistoGeneTroisCutNoOverF("higgsNN1_"+prefixe,"dipho_mgg",90,150,60,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN2_"+prefixe,"dipho_mgg",80,150,35,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN3_"+prefixe,"dipho_mgg",90,150,20,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN4_"+prefixe,"dipho_mgg",90,150,15,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN5_"+prefixe,"dipho_mgg",80,150,14,allTheCutsNN,myFile);	
*/

        TString cutEB ="&&pholead_isEB==1&&photrail_isEB==1";
        TString cutEE = "(pholead_isEE==1&&photrail_isEE==1)";
        TString cutMixte = "||((pholead_isEB==1&&photrail_isEE==1)||(pholead_isEE==1&&photrail_isEB==1))";
        TString cutNotEB = "&&(!(pholead_isEB==1&&photrail_isEB==1))";

        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat0","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat1","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9<0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat2","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat3","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9<0.94",myFile);
     
 
        if(theSelecType == 0){
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat0NN","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.7",myFile);
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat1NN","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.51",myFile);
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat2NN","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.6",myFile);
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat3NN","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.31",myFile);
        } 


        if(theSelecType == 1){
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat0NN","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.61",myFile);
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat1NN","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.37",myFile);
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat2NN","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.",myFile);
            doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat3NN","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.26",myFile);
        }



     /*
      doHistoGeneTroisCut("minr9","dipho_minR9",0,1.2,40,allTheCuts,myFile);
      doHistoGeneTroisCut("minr9TwoEB","dipho_minR9",0,1.2,40,allTheCuts+cutEB,myFile);
      doHistoGeneTroisCut("minr9TwoNotEB","dipho_minR9",0,1.2,40,allTheCuts+cutNotEB,myFile);

      doHistoGeneTroisCut("theCat","dipho_placeCats",0,1,2, allTheCuts,myFile);
      doHistoGeneTroisCut("theCatlowR9","dipho_placeCats",0,1,2, allTheCuts+"&&dipho_minR9<0.94",myFile);
      doHistoGeneTroisCut("theCathighR9","dipho_placeCats",0,1,2, allTheCuts+"&&dipho_minR9>0.94",myFile);
*/
myFile->Close();  

}
