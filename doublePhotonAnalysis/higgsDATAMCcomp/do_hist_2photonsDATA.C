#include "TString.h"
#include "TH1F.h"
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
//#include <iostream>

TString name = "data";

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

int printCutFlow(TString *cuts, int nbCuts, TString nameTag, TString localCut){

TH1F *h = new TH1F("h","",100,-3,3);
chain->Draw("pholead_eta>>h");
cout << "<" << nameTag << ">";
cout << "<cut0>" << h->Integral() << "</cut0>";
TString theCuts = "";
theCuts = localCut+"&&";
for (int i = 0 ; i < nbCuts; i++){
	theCuts+=cuts[i];
	chain->Draw("pholead_eta>>h",theCuts);
	cout << "<cut" << i+1 << ">" << h->Integral() << "</cut" << i+1 << ">";
	theCuts+="&&";
}
delete h;
cout << "</" << nameTag << ">";
}

int printAllCutFlow(TString *cuts, int nbCuts){
	cout << "<cutFlow_"+name+">" << endl;
	printCutFlow(cuts,nbCuts,"all","1");
	for (int i = 110; i < 160 ; i=i+10){
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


TH1F *th1fMaker(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts,  TString nomEnPlus, TString cutEnPlus ){
	TString theCut;
	/*if (part==1) {theCut = "pho_isEB==1"; nom+="_EB";}
	else if (part==2) {theCut = "pho_isEE==1"; nom+="_EE";} 
	else {theCut = "1"; nom+="_All";}*/
	theCut += Cuts+"&&"+cutEnPlus;
	//nom+="_"+nomEnPlus;
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

int doHistoGeneTroisCut(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts,  TFile *theFile){
	cout << "processing " << nom << endl;
	TH1F *data = (TH1F*) th1fMaker(nom,variable,Xmin,Xmax,bining,Cuts, "1", "1");
 	theFile->Write();
	delete data;
}

TH1F *th1fMakerNoOverF(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts,  TString nomEnPlus, TString cutEnPlus ){
	TString theCut;
	/*if (part==1) {theCut = "pho_isEB==1"; nom+="_EB";}
	 else if (part==2) {theCut = "pho_isEE==1"; nom+="_EE";} 
	 else {theCut = "1"; nom+="_All";}*/
	theCut += Cuts+"&&"+cutEnPlus;
	//nom+="_"+nomEnPlus;
	TH1F *h = new TH1F(nom,"",bining, Xmin, Xmax);
	TString toDraw = variable + ">>" + nom;
	chain->Draw(toDraw,theCut);
	return h;
}

int doHistoGeneTroisCutNoOverF(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts,  TFile *theFile){
	cout << "processing " << nom << endl;
	TH1F *data = (TH1F*) th1fMakerNoOverF(nom,variable,Xmin,Xmax,bining,Cuts, "1", "1");
 	theFile->Write();
	delete data;
}

TH1F *doNminus1(TString nom, TString variable, float Xmin, float Xmax, float bining, TString* cuts, int nbCut, int cutNb, TString nomEnPlus, TString cutEnPlus){
	TString theNminus1Cut= "1";
/*	if (part==1) {theNminus1Cut = "pho_isEB==1"; nom+="_EB";}
	else if (part==2) {theNminus1Cut = "pho_isEE==1"; nom+="_EE";} 
	else {theNminus1Cut = "1"; nom+="_All";}*/
	theNminus1Cut += "&&"+cutEnPlus;
	//nom+="_"+nomEnPlus;
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

int doHistoGeneTroisCutNminus1(TString nom, TString variable, float Xmin, float Xmax, float bining, TString* cuts, int nbCut, int cutNb, TFile *theFile){
        cout << "processing the plot " << nom << endl;
	TH1F *data = (TH1F*) doNminus1(nom,variable,Xmin,Xmax,bining, cuts, nbCut, cutNb,"1","1");
 	theFile->Write();

	delete data;
}

int do_hist_2photonsDATA(){
  TFile *myFile = new TFile("diphoFile_DATA_1.root","RECREATE");
   chain->Add("diphoton_part1.root");
  
  
//TString preSelection = "(abs(pholead_SCeta)<=2.5)&&(abs(photrail_SCeta)<=2.5)&&pholead_pt>30&&photrail_pt>30";
//TString etSelection = "pholead_pt>30&&photrail_pt>30";
int theSelecType=1;
TString prefixe;
int nbCuts = 6;
TString selection[6];
if (theSelecType==0){
prefixe = "CIC";
selection[0] = "(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&dipho_mgg>85";
selection[1] =  "(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))";
selection[2] =  "(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))";
selection[3] = "pholead_pt>40";
selection[4] = "photrail_pt>30";
//selection[5] = "pholead_pt>30 && (abs(pholead_eta)<=2.5)&&(!((abs(pholead_eta)>1.4442)&&(abs(pholead_eta)<1.566))) &&  pholead_TrackerIsodR03<(7.0+0.002*pholead_pt)&&photrail_EcalIsodR03<(8.2+0.012*pholead_pt)&&pholead_HcalIsodR03<(4.4+0.005*pholead_pt)&&((pholead_isEB==1&&pholead_sigieta<0.013)||(pholead_isEE==1&&pholead_sigieta<0.034))&&pholead_hoe<0.1&&photrail_pt>30 && (abs(photrail_eta)<=2.5)&&(!((abs(photrail_eta)>1.4442)&&(abs(photrail_eta)<1.566))) &&  photrail_TrackerIsodR03<(7.0+0.002*photrail_pt)&&photrail_EcalIsodR03<(8.2+0.012*photrail_pt)&&pholead_HcalIsodR03<(4.4+0.005*photrail_pt)&&((photrail_isEB==1&&photrail_sigieta<0.013)||(photrail_isEE==1&&photrail_sigieta<0.034))&&photrail_hoe<0.1&&event_nVertex>4&&event_nVertex<8";
selection[5] = "((((pholead_isEB==1&&pholead_r9>0.94)&&(pholead_CICcombIsoRho<5.38&&pholead_CICtrackIso<2.95&&pholead_CICworstCombRho<5.07&&pholead_hoe<0.099&&pholead_sigieta<0.01097&&pholead_CICdR>1.0))||((pholead_isEB==1&&pholead_r9<0.94)&&(pholead_CICcombIsoRho<2.2&&pholead_CICtrackIso<2.2&&pholead_CICworstCombRho<3.4&&pholead_hoe<0.062&&pholead_sigieta<0.0097&&pholead_r9>0.36&&pholead_CICdR>0.062))||((pholead_isEE==1&&pholead_r9>0.94)&&(pholead_CICcombIsoRho<1.77&&pholead_CICtrackIso<2.3&&pholead_CICworstCombRho<3.9&&pholead_hoe<0.065&&pholead_sigieta<0.028&&pholead_CICdR>1.0))||((pholead_isEE==1&&pholead_r9<=0.94)&&(pholead_CICcombIsoRho<1.29&&pholead_CICtrackIso<1.45&&pholead_CICworstCombRho<1.84&&pholead_hoe<0.048&&pholead_sigieta<0.027&&pholead_r9>0.32&&pholead_CICdR>1.0)))&&(((photrail_isEB==1&&photrail_r9>0.94)&&(photrail_CICcombIsoRho<5.38&&photrail_CICtrackIso<2.95&&photrail_CICworstCombRho<5.07&&photrail_hoe<0.099&&photrail_sigieta<0.01097&&photrail_CICdR>1.0))||((photrail_isEB==1&&photrail_r9<0.94)&&(photrail_CICcombIsoRho<2.2&&photrail_CICtrackIso<2.2&&photrail_CICworstCombRho<3.4&&photrail_hoe<0.062&&photrail_sigieta<0.0097&&photrail_r9>0.36&&photrail_CICdR>0.062))||((photrail_isEE==1&&photrail_r9>0.94)&&(photrail_CICcombIsoRho<1.77&&photrail_CICtrackIso<2.3&&photrail_CICworstCombRho<3.9&&photrail_hoe<0.065&&photrail_sigieta<0.028&&photrail_CICdR>1.0))||((photrail_isEE==1&&photrail_r9<0.94)&&(photrail_CICcombIsoRho<1.29&&photrail_CICtrackIso<1.45&&photrail_CICworstCombRho<1.84&&photrail_hoe<0.048&&photrail_sigieta<0.027&&photrail_r9>0.32&&photrail_CICdR>1.0))))";
// select CIC selection[5] = "((pholead_isEB==1&&(((pholead_r9>0.94)&&(pholead_CICcombIsoRho<3.8&&pholead_CICtrackIso<3.5&&pholead_CICworstCombRho<11.7&&pholead_hoe<0.082&&pholead_sigieta<0.0105&&pholead_CICdR>1.0))||((pholead_r9<0.94)&&(pholead_CICcombIsoRho<2.2&&pholead_CICtrackIso<2.2&&pholead_CICworstCombRho<3.4&&pholead_hoe<0.062&&pholead_sigieta<0.0097&&pholead_r9>0.36&&pholead_CICdR>0.062&&pholead_CICdR>1.0))))||((pholead_isEE==1&&(((pholead_r9>0.94)&&(pholead_CICcombIsoRho<1.77&&pholead_CICtrackIso<2.3&&pholead_CICworstCombRho<3.9&&pholead_hoe<0.065&&pholead_sigieta<0.028&&pholead_CICdR>1.0))||((pholead_r9<=0.94)&&(pholead_CICcombIsoRho<1.29&&pholead_CICtrackIso<1.45&&pholead_CICworstCombRho<1.84&&pholead_hoe<0.048&&pholead_sigieta<0.027&&pholead_r9>0.32&&pholead_CICdR>1.0))))))&&((photrail_isEB==1&&(((photrail_r9>0.94)&&(photrail_CICcombIsoRho<3.8&&photrail_CICtrackIso<3.5&&photrail_CICworstCombRho<11.7&&photrail_hoe<0.082&&photrail_sigieta<0.0105&&photrail_CICdR>1.0))||((photrail_r9<0.94)&&(photrail_CICcombIsoRho<2.2&&photrail_CICtrackIso<2.2&&photrail_CICworstCombRho<3.4&&photrail_hoe<0.062&&photrail_sigieta<0.0097&&photrail_r9>0.36&&photrail_CICdR>0.062))))||((photrail_isEE==1&&(((photrail_r9>0.94)&&(photrail_CICcombIsoRho<1.77&&photrail_CICtrackIso<2.3&&photrail_CICworstCombRho<3.9&&photrail_hoe<0.065&&photrail_sigieta<0.028&&photrail_CICdR>1.0))||((photrail_r9<0.94)&&(photrail_CICcombIsoRho<1.29&&photrail_CICtrackIso<1.45&&photrail_CICworstCombRho<1.84&&photrail_hoe<0.048&&photrail_sigieta<0.027&&photrail_r9>0.32&&photrail_CICdR>1.0))))))";
//selection[5] =  "((pholead_isEB==1&&(((pholead_r9>0.94)&&(pholead_CICcombIsoRho<3.33281&&pholead_CICtrackIso<8.55045&&pholead_CICworstCombRho<3.69087&&pholead_hoe<0.0970634&&pholead_sigieta<0.0106105&&pholead_CICdR>1.0))||((pholead_r9<0.94)&&(pholead_CICcombIsoRho<2.19529&&pholead_CICtrackIso<2.18377&&pholead_CICworstCombRho<2.42234&&pholead_hoe<0.0333386&&pholead_sigieta<0.0101019&&pholead_r9>0.283491&&pholead_CICdR>0.062&&pholead_CICdR>1.0))))||((pholead_isEE==1&&(((pholead_r9>0.94)&&(pholead_CICcombIsoRho<3.77659&&pholead_CICtrackIso<1.94267&&pholead_CICworstCombRho<6.41493&&pholead_hoe<0.062406&&pholead_sigieta<0.0290648&&pholead_CICdR>1.0))||((pholead_r9<=0.94)&&(pholead_CICcombIsoRho<1.63838&&pholead_CICtrackIso<5.48653&&pholead_CICworstCombRho<1.78649&&pholead_hoe<0.0741057&&pholead_sigieta<0.0270351&&pholead_r9>0.171847))))))&&((photrail_isEB==1&&(((photrail_r9>0.94)&&(photrail_CICcombIsoRho<3.33281&&photrail_CICtrackIso<8.55045&&photrail_CICworstCombRho<3.69087&&photrail_hoe<0.0970634&&photrail_sigieta<0.0106105&&photrail_CICdR>1.0))||((photrail_r9<0.94)&&(photrail_CICcombIsoRho<2.19529&&photrail_CICtrackIso<2.18377&&photrail_CICworstCombRho<2.42234&&photrail_hoe<0.0333386&&photrail_sigieta<0.0101019&&photrail_r9>0.283491&&photrail_CICdR>0.062))))||((photrail_isEE==1&&(((photrail_r9>0.94)&&(photrail_CICcombIsoRho<3.77659&&photrail_CICtrackIso<1.94267&&photrail_CICworstCombRho<6.41493&&photrail_hoe<0.062406&&photrail_sigieta<0.0290648&&photrail_CICdR>1.0))||((photrail_r9<0.94)&&(photrail_CICcombIsoRho<1.63838&&photrail_CICtrackIso<5.48653&&photrail_CICworstCombRho<1.78649&&photrail_hoe<0.0741057&&photrail_sigieta<0.0270351&&photrail_r9>0.171847&&photrail_CICdR>1.0))))))";
//selection[5] =  "pholead_isPassingCICRho==1&&photrail_isPassingCICRho==1";
}
/*
     TString cutEB ="&&pholead_isEB==1&&photrail_isEB==1";
        TString cutEE = "(pholead_isEE==1&&photrail_isEE==1)";
       TString cutMixte = "||((pholead_isEB==1&&photrail_isEE==1)||(pholead_isEE==1&&photrail_isEB==1))";
        TString cutNotEB = "&&(!(pholead_isEB==1&&photrail_isEB==1))";


TString CICvars[14] = {
"(((pholead_CICcombIso*50/pholead_pt)-0.17*event_rho)*pholead_pt/50)","(((photrail_CICcombIso*50/photrail_pt)-0.17*event_rho)*photrail_pt/50)",
"(((pholead_CICworstComb*50/pholead_pt)-0.52*event_rho)*pholead_pt/50)","(((photrail_CICworstComb*50/photrail_pt)-0.52*event_rho)*photrail_pt/50)",
"pholead_CICtrackIso","photrail_CICtrackIso",
"pholead_sigieta","photrail_sigieta",
"pholead_hoe","photrail_hoe",
"pholead_r9","photrail_r9",
"pholead_CICdR","photrail_CICdR"
};

TString CICnom[14] = {
"pholead_CICcombIso","photrail_CICcombIso",
"pholead_CICworstComb","photrail_CICworstComb",
"pholead_CICtrackIso","photrail_CICtrackIso",
"pholead_sigieta","photrail_sigieta",
"pholead_hoe","photrail_hoe",
"pholead_r9","photrail_r9",
"pholead_CICdR","photrail_CICdR"
};


float vcicST[7][4] = {
  {3.8,     2.2,     1.77,   1.29},    // rel combIso (good vtx)                                
  {11.7,    3.4,     3.9,    1.84},    // rel combIso (bad vtx)                                                   
  {3.5,     2.2,     2.3,    1.45},    // trkIso (good vtx)                                                           
  {0.0105,  0.0097,  0.028,  0.027},  // sigma_ie_ie                                                                  
  {0.082,   0.062,   0.065,  0.048},  // H/E       
  {0.94,    0.36,    0.94,   0.32},   // R9
  {1.,      0.062,   0.97,   0.97},    // dR to trk
};
:
TString NM1Cuts[7];
for (int k=0 ; k < 7 ; k++){
TString theCutCIC = "(";
for (int i =0 ; i < 4 ; i++){
	TString theCat = IntToString(i);
	theCutCIC+= "(pholead_Cat=="+theCat;
	for (int j = 0 ; j < 7 ; j++){
		if (j==k) continue; 
		TString theCut = FloatToString(vcicST[j][i]);
		TString comparateur = "<";
		if (j>4) comparateur = ">";
		theCutCIC+= "&&("+CICvars[2*j] + comparateur + theCut+")";
		cout << endl;
	}
	theCutCIC+=")||";
}
theCutCIC+="0)&&(";
for (int i =0 ; i < 4 ; i++){
        TString theCat = IntToString(i);
        theCutCIC+= "(photrail_Cat=="+theCat;
        for (int j = 0 ; j < 7 ; j++){
                if (j==k) continue; 
                TString theCut = FloatToString(vcicST[j][i]);
                TString comparateur = "<";
                if (j>4) comparateur = ">";
                theCutCIC+= "&&("+CICvars[2*j+1] + comparateur + theCut+")";
                cout << endl;
        }   
        theCutCIC+=")||";
}
theCutCIC+="0)";
cout << theCutCIC << endl;
NM1Cuts[k] = "(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&dipho_mgg>85&&(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))&&(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))&&pholead_pt>40&&photrail_pt>30&&"+theCutCIC;
}

TString theWholeCIC = "(";
for (int i =0 ; i < 4 ; i++){
        TString theCat = IntToString(i);
        theWholeCIC+= "(pholead_Cat=="+theCat;
        for (int j = 0 ; j < 7 ; j++){
                TString theCut = FloatToString(vcicST[j][i]);
                TString comparateur = "<";
                if (j>4) comparateur = ">";
                theWholeCIC+= "&&("+CICvars[2*j] + comparateur + theCut+")";
                cout << endl;
        }   
        theWholeCIC+=")||";
}
theWholeCIC+="0)&&(";
for (int i =0 ; i < 4 ; i++){
        TString theCat = IntToString(i);
        theWholeCIC+= "(photrail_Cat=="+theCat;
        for (int j = 0 ; j < 7 ; j++){
                TString theCut = FloatToString(vcicST[j][i]);
                TString comparateur = "<";
                if (j>4) comparateur = ">";
                theWholeCIC+= "&&("+CICvars[2*j+1] + comparateur + theCut+")";
                cout << endl;
        }
        theWholeCIC+=")||";
}
theWholeCIC+="0)";

cout << "chole " << theWholeCIC << endl;
TString totCut = "(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&dipho_mgg>85&&(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))&&(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))&&pholead_pt>40&&photrail_pt>30&&"+ theWholeCIC;
cout << "youhou the tot cut " << totCut << endl; 
doHistoGeneTroisCutNoOverF("higgs2_rhoCor","dipho_mgg",80,150,35,totCut,myFile);
return;
float catAll[21] = {
-2,6,32, //comb Iso
-5,15,80, //worst comb iso
0,3,30, // track iso 
0,0.1,80, //sig ieta
0,1.13,34, // HoE
0.3,1,35, // R9
0,0.07,35
};

float catCat0[21] = {
-2,6,32, //comb Iso
-5,15,80, //worst comb iso
0,3,30, // track iso 
0.005,0.016,33, //sig ieta
0,1.13,34, // HoE
0.94,1,30, // R9
0,0.07,35
};


float catCat1[21] = {
-2,6,32, //comb Iso
-5,15,80, //worst comb iso
0,3,30, // track iso 
0.005,0.016,33, //sig ieta
0,1.13,34, // HoE
0.3,0.94,32, // R9
0,0.07,35
};

float catCat2[21] = {
-2,6,32, //comb Iso
-5,15,80, //worst comb iso
0,3,30, // track iso 
0.018,0.038,40, //sig ieta
0,1.13,34, // HoE
0.94,1,30, // R9
0,0.07,35
};


float catCat3[21] = {
-2,6,32, //comb Iso
-5,15,80, //worst comb iso
0,3,30, // track iso 
0.018,0.038,40, //sig ieta
0,1.13,34, // HoE
0.3,0.94,32, // R9
0,0.07,35
};



for (int i = 0 ; i < 7 ; i++){
	doHistoGeneTroisCut(CICnom[2*i], CICvars[2*i],catAll[3*i],catAll[3*i+1], catAll[3*i+2],  NM1Cuts[i], myFile);
	doHistoGeneTroisCut(CICnom[2*i]+"_cat0", CICvars[2*i],catCat0[3*i],catCat0[3*i+1], catCat0[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==0", myFile);
	doHistoGeneTroisCut(CICnom[2*i]+"_cat1", CICvars[2*i],catCat1[3*i],catCat1[3*i+1], catCat1[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==1", myFile);
	doHistoGeneTroisCut(CICnom[2*i]+"_cat2", CICvars[2*i],catCat2[3*i],catCat2[3*i+1], catCat2[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==2", myFile);
	doHistoGeneTroisCut(CICnom[2*i]+"_cat3", CICvars[2*i],catCat3[3*i],catCat3[3*i+1], catCat3[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==3", myFile);


	doHistoGeneTroisCut(CICnom[2*i+1], CICvars[2*i+1],catAll[3*i],catAll[3*i+1], catAll[3*i+2],  NM1Cuts[i], myFile);
	doHistoGeneTroisCut(CICnom[2*i+1]+"_cat0", CICvars[2*i+1],catCat0[3*i],catCat0[3*i+1], catCat0[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==0", myFile);
	doHistoGeneTroisCut(CICnom[2*i+1]+"_cat1", CICvars[2*i+1],catCat1[3*i],catCat1[3*i+1], catCat1[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==1", myFile);
	doHistoGeneTroisCut(CICnom[2*i+1]+"_cat2", CICvars[2*i+1],catCat2[3*i],catCat2[3*i+1], catCat2[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==2", myFile);
	doHistoGeneTroisCut(CICnom[2*i+1]+"_cat3", CICvars[2*i+1],catCat3[3*i],catCat3[3*i+1], catCat3[3*i+2],  NM1Cuts[i]+"&&pholead_Cat==3", myFile);
}
return;*/

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
		0,1,50,  // lamdba ratio
		0,1,50,
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
		//	doHistoGeneTroisCut(Form("higgsTest_Cut%i",i),"dipho_mgg",90,150,60,allTheCuts,myFile);
//			doHistoGeneTroisCut(nom,varsToPlot[j],bining[j*3],bining[j*3+1],bining[j*3+2],allTheCuts,myFile);

			cout << "apres" << endl;
		}
		
	}
	//doHistoGeneTroisCut("vertex","event_nVertex",0,20,20,allTheCuts,myFile);
	
//	doHistoGeneTroisCutNoOverF("higgs1_"+prefixe,"dipho_mgg",70,150,80,allTheCuts,myFile);
	doHistoGeneTroisCutNoOverF("higgs2_"+prefixe,"dipho_mgg",80,150,35,allTheCuts,myFile);
//	doHistoGeneTroisCutNoOverF("higgs3_"+prefixe,"dipho_mgg",90,150,20,allTheCuts,myFile);
//	doHistoGeneTroisCutNoOverF("higgs4_"+prefixe,"dipho_mgg",90,150,15,allTheCuts,myFile);
//	doHistoGeneTroisCutNoOverF("higgs5_"+prefixe,"dipho_mgg",80,150,14,allTheCuts,myFile);

	
/*	TString allTheCutsNN = allTheCuts + "&&dipho_minNNshape>0.7";
	doHistoGeneTroisCutNoOverF("higgsNN1_"+prefixe,"dipho_mgg",90,150,60,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN2_"+prefixe,"dipho_mgg",80,150,35,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN3_"+prefixe,"dipho_mgg",90,150,20,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN4_"+prefixe,"dipho_mgg",90,150,15,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsNN5_"+prefixe,"dipho_mgg",80,150,14,allTheCutsNN,myFile);	
	
	
	TString allTheCutsNN = allTheCuts + "&&dipho_minNNshape>0.7&&(pholead_isEB==1&&photrail_isEB==1)&&dipho_minR9>0.94";
	doHistoGeneTroisCutNoOverF("higgsCat0NN1_"+prefixe,"dipho_mgg",90,150,60,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat0NN2_"+prefixe,"dipho_mgg",80,150,35,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat0NN3_"+prefixe,"dipho_mgg",90,150,20,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat0NN4_"+prefixe,"dipho_mgg",90,150,15,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat0NN5_"+prefixe,"dipho_mgg",80,150,14,allTheCutsNN,myFile);	

	TString allTheCutsNN = allTheCuts + "&&dipho_minNNshape>0.7&&(pholead_isEB==1&&photrail_isEB==1)&&dipho_minR9<0.94";
	doHistoGeneTroisCutNoOverF("higgsCat1NN1_"+prefixe,"dipho_mgg",90,150,60,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat1NN2_"+prefixe,"dipho_mgg",80,150,35,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat1NN3_"+prefixe,"dipho_mgg",90,150,20,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat1NN4_"+prefixe,"dipho_mgg",90,150,15,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat1NN5_"+prefixe,"dipho_mgg",80,150,14,allTheCutsNN,myFile);	


	TString allTheCutsNN = allTheCuts + "&&dipho_minNNshape>0.7&&!(pholead_isEB==1&&photrail_isEB==1)&&dipho_minR9>0.94";
	doHistoGeneTroisCutNoOverF("higgsCat2NN1_"+prefixe,"dipho_mgg",90,150,60,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat2NN2_"+prefixe,"dipho_mgg",80,150,35,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat2NN3_"+prefixe,"dipho_mgg",90,150,20,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat2NN4_"+prefixe,"dipho_mgg",90,150,15,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat2NN5_"+prefixe,"dipho_mgg",80,150,14,allTheCutsNN,myFile);	


	TString allTheCutsNN = allTheCuts + "&&dipho_minNNshape>0.7&&!(pholead_isEB==1&&photrail_isEB==1)&&dipho_minR9<0.94";
	doHistoGeneTroisCutNoOverF("higgsCat3NN1_"+prefixe,"dipho_mgg",90,150,60,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat3NN2_"+prefixe,"dipho_mgg",80,150,35,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat3NN3_"+prefixe,"dipho_mgg",90,150,20,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat3NN4_"+prefixe,"dipho_mgg",90,150,15,allTheCutsNN,myFile);
	doHistoGeneTroisCutNoOverF("higgsCat3NN5_"+prefixe,"dipho_mgg",80,150,14,allTheCutsNN,myFile);	



*/

     TString cutEB ="&&pholead_isEB==1&&photrail_isEB==1";
        TString cutEE = "(pholead_isEE==1&&photrail_isEE==1)";
       TString cutMixte = "||((pholead_isEB==1&&photrail_isEE==1)||(pholead_isEE==1&&photrail_isEB==1))";
        TString cutNotEB = "&&(!(pholead_isEB==1&&photrail_isEB==1))";

  /*      doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat0","dipho_mgg",100,150,25,allTheCuts+cutEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat1","dipho_mgg",100,150,25,allTheCuts+cutEB+"&&dipho_minR9<0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat2","dipho_mgg",100,150,25,allTheCuts+cutNotEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"_cat3","dipho_mgg",100,150,25,allTheCuts+cutNotEB+"&&dipho_minR9<0.94",myFile);
*/
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"L_cat0","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"L_cat1","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9<0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"L_cat2","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"L_cat3","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9<0.94",myFile);
/*	
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"B_cat0","dipho_mgg",100,150,200,allTheCuts+cutEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"B_cat1","dipho_mgg",100,150,200,allTheCuts+cutEB+"&&dipho_minR9<0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"B_cat2","dipho_mgg",100,150,200,allTheCuts+cutNotEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("higgs2_"+prefixe+"B_cat3","dipho_mgg",100,150,200,allTheCuts+cutNotEB+"&&dipho_minR9<0.94",myFile);


        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"_cat0","dipho_mgg",100,150,25,allTheCuts+cutEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.65",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"_cat1","dipho_mgg",100,150,25,allTheCuts+cutEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.30",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"_cat2","dipho_mgg",100,150,25,allTheCuts+cutNotEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.15",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"_cat3","dipho_mgg",100,150,25,allTheCuts+cutNotEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.12",myFile);

        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"L_cat0","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.65",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"L_cat1","dipho_mgg",80,150,35,allTheCuts+cutEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.30",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"L_cat2","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.15",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"L_cat3","dipho_mgg",80,150,35,allTheCuts+cutNotEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.12",myFile);

        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"B_cat0","dipho_mgg",100,150,200,allTheCuts+cutEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.65",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"B_cat1","dipho_mgg",100,150,200,allTheCuts+cutEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.30",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"B_cat2","dipho_mgg",100,150,200,allTheCuts+cutNotEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.15",myFile);
        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe+"B_cat3","dipho_mgg",100,150,200,allTheCuts+cutNotEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.12",myFile);



        doHistoGeneTroisCutNoOverF("higgs2NN_"+prefixe,"dipho_mgg",80,150,35,allTheCuts+"&&(1"+cutEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.65||1"+cutEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.30||1"+cutNotEB+"&&dipho_minR9>0.94&&dipho_minNNshape>0.15||1"+cutNotEB+"&&dipho_minR9<0.94&&dipho_minNNshape>0.12)",myFile);

        doHistoGeneTroisCutNoOverF("minNN_"+prefixe,"dipho_minNNshape",-0.2,1.2,40,allTheCuts,myFile);

        doHistoGeneTroisCutNoOverF("minNN_"+prefixe+"_cat0","dipho_minNNshape",-0.2,1.2,40,allTheCuts+cutEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("minNN_"+prefixe+"_cat1","dipho_minNNshape",-0.2,1.2,40,allTheCuts+cutEB+"&&dipho_minR9<0.94",myFile);
        doHistoGeneTroisCutNoOverF("minNN_"+prefixe+"_cat2","dipho_minNNshape",-0.2,1.2,40,allTheCuts+cutNotEB+"&&dipho_minR9>0.94",myFile);
        doHistoGeneTroisCutNoOverF("minNN_"+prefixe+"_cat3","dipho_minNNshape",-0.2,1.2,40,allTheCuts+cutNotEB+"&&dipho_minR9<0.94",myFile);


/*      doHistoGeneTroisCut("minr9","dipho_minR9",0,1.2,40,allTheCuts,myFile);
      doHistoGeneTroisCut("minr9TwoEB","dipho_minR9",0,1.2,40,allTheCuts+cutEB,myFile);
      doHistoGeneTroisCut("minr9TwoNotEB","dipho_minR9",0,1.2,40,allTheCuts+cutNotEB,myFile);

      doHistoGeneTroisCut("theCat","dipho_placeCats",0,1,2, allTheCuts,myFile);
      doHistoGeneTroisCut("theCatlowR9","dipho_placeCats",0,1,2, allTheCuts+"&&dipho_minR9<0.94",myFile);
      doHistoGeneTroisCut("theCathighR9","dipho_placeCats",0,1,2, allTheCuts+"&&dipho_minR9>0.94",myFile);


*/

myFile->Close();  




}
