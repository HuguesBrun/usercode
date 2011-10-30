#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

float signifToUse = 1140;

TFile *myFile = new TFile("theSignifile2.root","RECREATE");
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
TString FloatToString(float number){
  ostringstream oss;
  oss << number;
  return oss.str();
}


// Higgs signal
TChain *HiggsSignalGluGlu = new TChain("diPhotons");
TChain *HiggsSignalHiggstralung  = new TChain("diPhotons");
TChain *HiggsTTfusion = new TChain("diPhotons");
TChain *HiggsSignalVBF  = new TChain("diPhotons");


//reducible bg
TChain *QCD40doubleEMenriched  = new TChain("diPhotons");
TChain *QCD30to40doubleEMenriched = new TChain("diPhotons");
TChain *GJet =  new TChain("diPhotons");
TChain *DYtoEE = new TChain("diPhotons");

//irreducible bg
TChain *diPhotonBox10to25 = new TChain("diPhotons");
TChain *diPhotonBox25to250 = new TChain("diPhotons");
TChain *diPhotonBox250toInf = new TChain("diPhotons");
TChain *diPhotonBorn10to25 = new TChain("diPhotons");
TChain *diPhotonBorn25to250 = new TChain("diPhotons");
TChain *diPhotonBorn250toInf = new TChain("diPhotons");

vector <float> calcNumber(TString cuts, float intLumi){
	TString cutDiphoton = "((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))";
	TString cutGJet = "((((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1))&&(!((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))||((!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	TString cutDiJet = "(!((pholead_isPromptGenPho==1)||(pholead_isFromQuarkGen==1)))&&(!(((photrail_isPromptGenPho==1)||(photrail_isFromQuarkGen==1))))";
	
	vector <float> result;
	
	float nbSignal = 
	HiggsSignalGluGlu->GetEntries(cuts)*0.038617/109989
	+ HiggsSignalHiggstralung->GetEntries(cuts)*0.002482/102290
	+ HiggsTTfusion->GetEntries(cuts)*0.000236/22000
	+ HiggsSignalVBF->GetEntries(cuts)*0.000236/22000;
	result.push_back(nbSignal*intLumi);
	float nbDipho = 
	GJet->GetEntries(cuts+"&&"+cutDiphoton)*501.15/6757937
	// diphoton box
	+ diPhotonBox10to25->GetEntries(cuts+"&&"+cutDiphoton)*358.20/528400
	+ diPhotonBox25to250->GetEntries(cuts+"&&"+cutDiphoton)*12.37/518288
	+ diPhotonBox250toInf->GetEntries(cuts+"&&"+cutDiphoton)*0.000208/515028
	// diphoton born
	+ diPhotonBorn10to25->GetEntries(cuts+"&&"+cutDiphoton)*236.40/505456
	+ diPhotonBorn25to250->GetEntries(cuts+"&&"+cutDiphoton)*22.37/532864
	+ diPhotonBorn250toInf->GetEntries(cuts+"&&"+cutDiphoton)*0.008072/52550
	// Di Jet 
	+ QCD30to40doubleEMenriched->GetEntries(cuts+"&&"+cutDiphoton)*10868.00/6074670
	+ QCD40doubleEMenriched->GetEntries(cuts+"&&"+cutDiphoton)*43571.00/39387002;
    result.push_back(nbDipho*intLumi*1.3);
   	cout << nbDipho*intLumi*1.3 << endl; 
	
	float nbGJet = GJet->GetEntries(cuts+"&&"+cutGJet)*501.15/6757937
	+ QCD30to40doubleEMenriched->GetEntries(cuts+"&&"+cutGJet)*10868.00/6074670
	+ QCD40doubleEMenriched->GetEntries(cuts+"&&"+cutGJet)*43571.00/39387002;
	result.push_back(nbGJet*intLumi*1.3);

	float nbDiJet = QCD40doubleEMenriched->GetEntries(cuts+"&&"+cutDiJet)*43571.00/39387002
	+ QCD30to40doubleEMenriched->GetEntries(cuts+"&&"+cutDiJet)*10868.00/6074670;
	result.push_back(nbDiJet*intLumi);
	
	float nbDY = DYtoEE->GetEntries(cuts+"&&"+cutDiJet)*1300.00/5312220;
	result.push_back(nbDY*intLumi*1.15);
	cout << nbDY*intLumi*1.15 << endl;
		
	return result;
	
}


TString giveMeTheCut(TMVA::MethodCuts* cutEBiR9, float effEBiR9){
//TString theLocalCut = "(dipho_HLT_bit6>-1||dipho_HLT_bit16>-1)&&(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&pholead_HasPixSeed==0&&photrail_HasPixSeed==0&&pholead_TrackerIso<2&&photrail_TrackerIso<2&&pholead_EcalIso<4.2&&photrail_EcalIso<4.2&&pholead_HcalIso<2.2&&photrail_HcalIso<2.2&&pholead_pt>40&&photrail_pt>30&&(((pholead_isEB==1)&&(pholead_sigieta<0.01))||((pholead_isEE==1)&&(pholead_sigieta<0.03)))&&(((photrail_isEB==1)&&(photrail_sigieta<0.01))||((photrail_isEE==1)&&(photrail_sigieta<0.03)))&&pholead_hoe<0.05&&photrail_hoe<0.05&&(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))&&(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))&&1";
std::vector<Double_t> cutsMinEBiR9;
std::vector<Double_t> cutsMaxEBiR9;
cutEBiR9->GetCuts( effEBiR9, cutsMinEBiR9, cutsMaxEBiR9 );


TString theLocalCutEBiR9 = "";
theLocalCutEBiR9 += "&&pholead_NNshapeOutput>"+FloatToString(cutsMinEBiR9[0]);


TString theLocalCutEBiR9Trail = "";
theLocalCutEBiR9Trail += "&&photrail_NNshapeOutput>"+FloatToString(cutsMinEBiR9[0]);




theLocalCut = theLocalCutEBiR9 + theLocalCutEBiR9Trail;

return theLocalCut; 
}


do_signif4CatEGM006NNbis(){
		using namespace RooFit;	
	HiggsSignalGluGlu->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_GluGluToHToGGM115.root");
	HiggsSignalHiggstralung->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_WHZHHToGGM115.root");
	HiggsTTfusion->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_TTHHToGGM115.root");
	HiggsSignalVBF->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_VBFHToGGM115.root");
	QCD40doubleEMenriched->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_QCDPt40doubleEMEnriched.root");
	QCD30to40doubleEMenriched->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_QCDPt30to40doubleEMEnriched.root");
	GJet->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_GJetPt20.root");
	diPhotonBox10to25->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_diphoBox10to25.root");
	diPhotonBox25to250->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_diphoBox25to250.root");
	diPhotonBox250toInf->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_diphoBox250toInf.root");
	diPhotonBorn10to25->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_diphoBorn10to25.root");
	diPhotonBorn25to250->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_diphoBorn25to250.root");
	diPhotonBorn250toInf->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_diphoBorn250toInf.root");
	DYtoEE->Add("/sps/cms/hbrun/DiPhotons42X/diPhotonMC/theDiphoMiniTrees/miniTreeEnriched/diphotonsEnriched_DY.root");

	
	TString preSelected[14] = {"pholead_TrackerIsodR03<(7.0+0.002*pholead_pt)",
		"photrail_TrackerIsodR03<(7.0+0.002*photrail_pt)",
		"pholead_EcalIsodR03<(8.2+0.012*pholead_pt)",
		"photrail_EcalIsodR03<(8.2+0.012*photrail_pt)",
		"pholead_HcalIsodR03<(4.4+0.005*pholead_pt)",
		"photrail_HcalIsodR03<(4.4+0.005*photrail_pt)",
		"((pholead_isEB==1&&pholead_sigieta<0.013)||(pholead_isEE==1&&pholead_sigieta<0.034))",
		"((photrail_isEB==1&&photrail_sigieta<0.013)||(photrail_isEE==1&&photrail_sigieta<0.034))",
		"pholead_hoe<0.1",
		"photrail_hoe<0.1",
		"pholead_pt>40",
		"photrail_pt>30",
		"pholead_HasPixSeed==0",
		"photrail_HasPixSeed==0"};
	
	TString preSelectedCuts = "";
	for (int i = 0 ; i < 12 ; i++){
			preSelectedCuts += preSelected[i] + "&&";
		
	}
        cout << preSelectedCuts << endl;
	preSelectedCuts += "1";
	TString selection[17];
	selection[0] = "(dipho_HLT_bit6>-1||dipho_HLT_bit16>-1)&&(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)";
	selection[1] =  "pholead_HasPixSeed==0";
	selection[2] =  "photrail_HasPixSeed==0";
	selection[3] =  "pholead_TrackerIso<2";
	selection[4] =  "photrail_TrackerIso<2";
	selection[5] =  "pholead_EcalIso<4.2";
	selection[6] =  "photrail_EcalIso<4.2";
	selection[7] =  "pholead_HcalIso<2.2";
	selection[8] =  "photrail_HcalIso<2.2";
	selection[9] =  "pholead_pt>40";
	selection[10] =         "photrail_pt>30";//30
	selection[11] =         "(((pholead_isEB==1)&&(pholead_sigieta<0.01))||((pholead_isEE==1)&&(pholead_sigieta<0.03)))";
	selection[12] =         "(((photrail_isEB==1)&&(photrail_sigieta<0.01))||((photrail_isEE==1)&&(photrail_sigieta<0.03)))";        
        selection[13] =         "pholead_hoe<0.05";
	selection[14] =         "photrail_hoe<0.05";
	selection[15] =         "(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))";
	selection[16] =         "(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))";
	TString selectionCut = "";
	for (int i = 0 ; i < 17 ; i++){
			selectionCut += selection[i] + "&&";
	}
	selectionCut += "1";
	TString selectionCut = "pholead_pt>40&&photrail_pt>30&&(!((abs(pholead_etaSC)>1.4442)&&(abs(pholead_etaSC)<1.566)))&&(!((abs(photrail_etaSC)>1.4442)&&(abs(photrail_etaSC)<1.566)))&&(abs(pholead_etaSC)<=2.5)&&(abs(photrail_etaSC)<=2.5)&&pholead_isPassingCICRho==1&&photrail_isPassingCICRho==1";

	cout << selectionCut << endl;
	TString cutEB ="&&pholead_isEB==1&&photrail_isEB==1";
	TString cutEE = "(pholead_isEE==1&&photrail_isEE==1)";
	TString cutMixte = "||((pholead_isEB==1&&photrail_isEE==1)||(pholead_isEE==1&&photrail_isEB==1))";
        TString cutNotEB = "&&(!(pholead_isEB==1&&photrail_isEB==1))";
//	giveNumbers(preSelectedCuts+"&&dipho_mgg>=90&&dipho_mgg<=150&&pholead_isEB==1&&photrail_isEB==1", 1000);
	cout << " all EGM loose 006" << endl;
//	giveNumbers(selectionCut+"&&dipho_mgg>=90&&dipho_mgg<=150", 1000);
//   	float theSigni = calcSignificance(selectionCut+"&&dipho_mgg>=113&&dipho_mgg<=117", 1000);
//	cout << "signi EGM LOOSE 006 " << theSigni << endl;


	// create the TMVA reader 

        TMVA::Reader *reader0 = new TMVA::Reader( "!Color:!Silent" );   
        float pho_IsoHollowTrkCone, pho_IsoEcalRechit, pho_IsoHcalRechit, pho_hoe, pho_sigmaIetaIeta, pho_theNN;
         reader0->AddVariable("photonsNN.pho_NNshapeOutput",&pho_theNN); 
	reader0->BookMVA("Cuts","/sps/cms/hbrun/theNN/NN42/higgsPrelec/weights_13/TMVAClassification_Cuts.weights.xml");
        TMVA::MethodCuts* mcuts0 = reader0->FindCutsMVA( "Cuts" ) ;       
  




	/*TString testCuts = giveMeTheCut(mcuts1,mcuts2,0.91,0.59);
	cout << testCuts << endl;
        cout << "theImproved signi = " << calcSignificance("dipho_mgg>=113&&dipho_mgg<=117&&"+testCuts, 1000) << endli;*/
 
	/*TString testCuts = giveMeTheCut(mcuts0, 0.95);
	cout << "les cuts = " << testCuts << endl;
        cout << "cat0 = " << calcSignificance("dipho_mgg>=100&&dipho_mgg<=150&&dipho_minR9>0.94&&"+testCuts+cutEB, 1000) << endl;
        cout << "cat1 = " << calcSignificance("dipho_mgg>=100&&dipho_mgg<=150&&dipho_minR9<=0.94&&"+testCuts+cutEB, 1000) << endl;
        cout << "cat2 = " << calcSignificance("dipho_mgg>=100&&dipho_mgg<=150&&dipho_minR9>0.94&&"+testCuts+cutNotEB, 1000) << endl;
        cout << "cat3 = " << calcSignificance("dipho_mgg>=100&&dipho_mgg<=150&&dipho_minR9<=0.94&&"+testCuts+cutNotEB, 1000) << endl;
        cout << "tot  = " << calcSignificance("dipho_mgg>=100&&dipho_mgg<=150&&"+testCuts, 1000) << endl;
*/
	
	
	float effCat0, effCat1, effCat2, effCat3;
	float effSig[100], effRedBg[100], signif[100], theEffCat[100], theNN[100];
	TString theCut;


        TString theCutEnPlus = "&&dipho_minR9>0.94&&dipho_mgg>=111&&dipho_mgg<=119"+cutEB;
	vector <float> beforeSelec = calcNumber(preSelectedCuts+theCutEnPlus, signifToUse);
	float numberSigPresel = beforeSelec[0];
	float numberBgPresel = beforeSelec[2]+beforeSelec[3]+beforeSelec[4];
	vector <float> afterSelec;
	for (int i=0 ; i < 100 ; i++){
		effCat0 = 1.0*i/100;
		theNN[i] = effCat0;
		theCut = "&&dipho_minNNshape>"+FloatToString(effCat0);
		cout << "coucou" << selectionCut+"&&"+theCut+theCutEnPlus << endl;
		afterSelec = calcNumber(selectionCut+theCut+theCutEnPlus, signifToUse);
		effSig[i] = afterSelec[0]/numberSigPresel;
		effRedBg[i] = 1-1.0*(afterSelec[2]+afterSelec[3]+afterSelec[4])/numberBgPresel;
		if ((afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4])>0 ) signif[i] = 1.0*afterSelec[0]/sqrt(afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4]);
		else signif[i] = 0;
		theEffCat[i] = effCat0;
		cout << "nb sig = " << afterSelec[0] << " nb bg = " << afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4] << endl;
		cout << "effSig=" << effSig[i] << " effBg=" << effRedBg[i] << " signif=" << signif[i] << endl; 
	}
	
	TGraph *graph1 = new TGraph(100, effSig, signif);
	graph1->SetName("sigCat0");
	graph1->Draw();
	graph1->Write();
	TGraph *graph2 = new TGraph(100, effSig, effRedBg);
	graph2->SetName("effCat0");
	graph2->Draw();
	graph2->Write();
	TGraph *graph3 = new TGraph(100, theEffCat, effSig);
	graph3->SetName("effSigCat0");
	graph3->Draw();
	graph3->Write();
	TGraph *graph4 = new TGraph(100,theEffCat, effRedBg);
	graph4->SetName("bgRegCat0");
	graph4->Draw();
	graph4->Write();
	TGraph *graph5 = new TGraph(100,theNN, signif);
	graph5->SetName("signifNNcat0");
	graph5->Draw();
	graph5->Write();
        TString theCutEnPlus = "&&dipho_minR9<0.94&&dipho_mgg>=111&&dipho_mgg<=119"+cutEB;
	vector <float> beforeSelec = calcNumber(preSelectedCuts+theCutEnPlus, signifToUse);
	float numberSigPresel = beforeSelec[0];
	float numberBgPresel = beforeSelec[2]+beforeSelec[3]+beforeSelec[4];
	vector <float> afterSelec;
	for (int i=0 ; i < 100 ; i++){
		effCat0 = 1.0*i/100;
		theNN[i] = effCat0;
		theCut = "&&dipho_minNNshape>"+FloatToString(effCat0);
		cout << "coucou" << selectionCut+"&&"+theCut+theCutEnPlus << endl;
		afterSelec = calcNumber(selectionCut+theCut+theCutEnPlus, signifToUse);
		effSig[i] = afterSelec[0]/numberSigPresel;
		effRedBg[i] = 1-1.0*(afterSelec[2]+afterSelec[3]+afterSelec[4])/numberBgPresel;
		if ((afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4])>0 ) signif[i] = 1.0*afterSelec[0]/sqrt(afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4]);
		else signif[i] = 0;
		theEffCat[i] = effCat0;
		cout << "effSig=" << effSig[i] << " effBg=" << effRedBg[i] << " signif=" << signif[i] << endl; 
		cout << "nb sig = " << afterSelec[0] << " nb bg = " << afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4] << endl;
	}
	
	TGraph *graph1 = new TGraph(100, effSig, signif);
	graph1->SetName("sigCat1");
	graph1->Draw();
	graph1->Write();
	TGraph *graph2 = new TGraph(100, effSig, effRedBg);
	graph2->SetName("effCat1");
	graph2->Draw();
	graph2->Write();
	TGraph *graph3 = new TGraph(100, theEffCat, effSig);
	graph3->SetName("effSigCat1");
	graph3->Draw();
	graph3->Write();
	TGraph *graph4 = new TGraph(100,theEffCat, effRedBg);
	graph4->SetName("bgRegCat1");
	graph4->Draw();
	graph4->Write();
	TGraph *graph5 = new TGraph(100,theNN, signif);
	graph5->SetName("signifNNcat1");
	graph5->Draw();
	graph5->Write();
        TString theCutEnPlus = "&&dipho_minR9>0.94&&dipho_mgg>=111&&dipho_mgg<=119"+cutNotEB;
	vector <float> beforeSelec = calcNumber(preSelectedCuts+theCutEnPlus, signifToUse);
	float numberSigPresel = beforeSelec[0];
	float numberBgPresel = beforeSelec[2]+beforeSelec[3]+beforeSelec[4];
	vector <float> afterSelec;
	for (int i=0 ; i < 100 ; i++){
		effCat0 = 1.0*i/100;
		theNN[i] = effCat0;
		theCut = "&&dipho_minNNshape>"+FloatToString(effCat0);
		cout << "coucou" << selectionCut+"&&"+theCut+theCutEnPlus << endl;
		afterSelec = calcNumber(selectionCut+theCut+theCutEnPlus, signifToUse);
		effSig[i] = afterSelec[0]/numberSigPresel;
		effRedBg[i] = 1-1.0*(afterSelec[2]+afterSelec[3]+afterSelec[4])/numberBgPresel;
		if ((afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4])>0 ) signif[i] = 1.0*afterSelec[0]/sqrt(afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4]);
		else signif[i] = 0;
		theEffCat[i] = effCat0;
		cout << "effSig=" << effSig[i] << " effBg=" << effRedBg[i] << " signif=" << signif[i] << endl; 
	}
	
	TGraph *graph1 = new TGraph(100, effSig, signif);
	graph1->SetName("sigCat2");
	graph1->Draw();
	graph1->Write();
	TGraph *graph2 = new TGraph(100, effSig, effRedBg);
	graph2->SetName("effCat2");
	graph2->Draw();
	graph2->Write();
	TGraph *graph3 = new TGraph(100, theEffCat, effSig);
	graph3->SetName("effSigCat2");
	graph3->Draw();
	graph3->Write();
	TGraph *graph4 = new TGraph(100,theEffCat, effRedBg);
	graph4->SetName("bgRegCat2");
	graph4->Draw();
	graph4->Write();
	TGraph *graph5 = new TGraph(100,theNN, signif);
	graph5->SetName("signifNNcat2");
	graph5->Draw();
	graph5->Write();

        TString theCutEnPlus = "&&dipho_minR9<0.94&&dipho_mgg>=111&&dipho_mgg<=119"+cutNotEB;
	vector <float> beforeSelec = calcNumber(preSelectedCuts+theCutEnPlus, signifToUse);
	float numberSigPresel = beforeSelec[0];
	float numberBgPresel = beforeSelec[2]+beforeSelec[3]+beforeSelec[4];
	vector <float> afterSelec;
	for (int i=0 ; i < 100 ; i++){
		effCat0 = 1.0*i/100;
		theNN[i] = effCat0;
		theCut = "&&dipho_minNNshape>"+FloatToString(effCat0);
		cout << "coucou" << selectionCut+"&&"+theCut+theCutEnPlus << endl;
		afterSelec = calcNumber(selectionCut+theCut+theCutEnPlus, signifToUse);
		effSig[i] = afterSelec[0]/numberSigPresel;
		effRedBg[i] = 1-1.0*(afterSelec[2]+afterSelec[3]+afterSelec[4])/numberBgPresel;
		if ((afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4])>0 ) signif[i] = 1.0*afterSelec[0]/sqrt(afterSelec[1]+afterSelec[2]+afterSelec[3]+afterSelec[4]);
		else signif[i] = 0;
		theEffCat[i] = effCat0;
		cout << "effSig=" << effSig[i] << " effBg=" << effRedBg[i] << " signif=" << signif[i] << endl; 
	}
	
	TGraph *graph1 = new TGraph(100, effSig, signif);
	graph1->SetName("sigCat3");
	graph1->Draw();
	graph1->Write();
	TGraph *graph2 = new TGraph(100, effSig, effRedBg);
	graph2->SetName("effCat3");
	graph2->Draw();
	graph2->Write();
	TGraph *graph3 = new TGraph(100, theEffCat, effSig);
	graph3->SetName("effSigCat3");
	graph3->Draw();
	graph3->Write();
	TGraph *graph4 = new TGraph(100,theEffCat, effRedBg);
	graph4->SetName("bgRegCat3");
	graph4->Draw();
	graph4->Write();
	TGraph *graph5 = new TGraph(100,theNN, signif);
	graph5->SetName("signifNNcat3");
	graph5->Draw();
	graph5->Write();
	myFile->Write();
	myFile->Close();
}
