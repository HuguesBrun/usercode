

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

do_theHistos(TString nom, TString variable, float Xmin, float Xmax, float bining, TString Cuts, int part, TString nomEnPlus, TString cutEnPlus ){
cout << "plotting " << variable << endl;
TH1F *h = th1fMaker(nom, variable, Xmin, Xmax, bining, Cuts, 0, nomEnPlus, cutEnPlus );
myFile->Write();
delete h;
TH1F *h = th1fMaker(nom, variable, Xmin, Xmax, bining, Cuts, 1, nomEnPlus, cutEnPlus );
myFile->Write();
delete h;
TH1F *h = th1fMaker(nom, variable, Xmin, Xmax, bining, Cuts, 2, nomEnPlus, cutEnPlus );
myFile->Write();
delete h;
}

do_scatPlot(TString nom, TString variable1, TString variable2, float Xmin, float Xmax, float binning, TString Cuts, int part, TString nomEnPlus, TString cutEnPlus) {

	TString localPart="1";
	if (part ==1) {localPart="pho_isEB==1"; nom +="_EB";}
	else if (part ==2) {localPart="pho_isEE==1"; nom += "_EE";}
	else nom += "_all";

	TString theCut = Cuts+"&&"+cutEnPlus+"&&"+localPart;
	TString theTotalNom = nom+"_"+nomEnPlus;
	TString toPlot = variable1+":"+variable2+">>"+theTotalNom;
	TH2F *scat = new TH2F(theTotalNom,"",binning, Xmin, Xmax, binning, Xmin, Xmax);
	chain->Draw(toPlot,theCut);
	myFile->Write();
	delete scat;

}

do_AllscatPlot(TString nom, TString variable1, TString variable2, float Xmin, float Xmax, float binning, TString Cuts, int part, TString nomEnPlus, TString cutEnPlus) {
		cout << "plotting " << variable1 << " versus " << variable2 << endl;	
		do_scatPlot(nom, variable1, variable2, Xmin, Xmax, binning, Cuts, 0, nomEnPlus, cutEnPlus);
		do_scatPlot(nom, variable1, variable2, Xmin, Xmax, binning, Cuts, 1, nomEnPlus, cutEnPlus);
		do_scatPlot(nom, variable1, variable2, Xmin, Xmax, binning, Cuts, 2, nomEnPlus, cutEnPlus);
}
TChain *chain = new TChain("photons");
TFile *myFile = new TFile("theIsoHisto.root","RECREATE");
do_IsoPlots(){
	chain->Add("MinTree3.root");
	TString theBaseCut="isAspike==0&&(abs(pho_SCeta)<=2.5)&&(!((abs(pho_SCeta)>1.4442)&&(abs(pho_SCeta)<1.566)))&&pho_hoe<0.15&&pho_et>18";
	
	TString theTriggerCutIsol = "pho_IsoEcalRechit03<(5+0.006*pho_et)&&pho_IsoHcalRechit03<(3+0.0025*pho_et)";
	TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(3+0.001*pho_et)&&pho_IsoEcalRechit03<(5+0.006*pho_et)&&pho_IsoHcalRechit03<(3+0.0025*pho_et)";
	

	do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"Isol",theTriggerCutIsol);
	do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"Isol",theTriggerCutIsol);
	
	do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"Isol",theTotalIdTriggerCutIsol);
	
	do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "Isol", theTriggerCutIsol);


///////////////////////////////////////////////////////////////
TString theTriggerCutIsol = "pho_IsoEcalRechit03<(5.5+0.006*pho_et)&&pho_IsoHcalRechit03<(3.5+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(3.5+0.001*pho_et)&&pho_IsoEcalRechit03<(5.5+0.006*pho_et)&&pho_IsoHcalRechit03<(3.5+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"LooseIsol",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"LooseIsol",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"LooseIsol",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "LooseIsol", theTriggerCutIsol);

///////////////////////////////////////////////////////////////////////

TString theTriggerCutIsol = "pho_IsoEcalRechit03<(6+0.006*pho_et)&&pho_IsoHcalRechit03<(4+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(4+0.001*pho_et)&&pho_IsoEcalRechit03<(6+0.006*pho_et)&&pho_IsoHcalRechit03<(4+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"LooserIsol",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"LooserIsol",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"LooserIsol",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "LooserIsol", theTriggerCutIsol);


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
TString theBaseCut="isAspike==0&&(abs(pho_SCeta)<=2.5)&&(!((abs(pho_SCeta)>1.4442)&&(abs(pho_SCeta)<1.566)))&&pho_hoe<0.15&&pho_et>22";

TString theTriggerCutIsol = "pho_IsoEcalRechit03<(5+0.006*pho_et)&&pho_IsoHcalRechit03<(3+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(3+0.001*pho_et)&&pho_IsoEcalRechit03<(5+0.006*pho_et)&&pho_IsoHcalRechit03<(3+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"IsolMidEt",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"IsolMidEt",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"IsolMidEt",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "IsolMidEt", theTriggerCutIsol);


///////////////////////////////////////////////////////////////
TString theTriggerCutIsol = "pho_IsoEcalRechit03<(5.5+0.006*pho_et)&&pho_IsoHcalRechit03<(3.5+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(3.5+0.001*pho_et)&&pho_IsoEcalRechit03<(5.5+0.006*pho_et)&&pho_IsoHcalRechit03<(3.5+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"LooseIsolMidEt",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"LooseIsolMidEt",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"LooseIsolMidEt",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "LooseIsolMidEt", theTriggerCutIsol);

///////////////////////////////////////////////////////////////////////

TString theTriggerCutIsol = "pho_IsoEcalRechit03<(6+0.006*pho_et)&&pho_IsoHcalRechit03<(4+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(4+0.001*pho_et)&&pho_IsoEcalRechit03<(6+0.006*pho_et)&&pho_IsoHcalRechit03<(4+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"LooserIsolMidEt",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"LooserIsolMidEt",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"LooserIsolMidEt",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "LooserIsolMidEt", theTriggerCutIsol);

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
TString theBaseCut="isAspike==0&&(abs(pho_SCeta)<=2.5)&&(!((abs(pho_SCeta)>1.4442)&&(abs(pho_SCeta)<1.566)))&&pho_hoe<0.15&&pho_et>26";

TString theTriggerCutIsol = "pho_IsoEcalRechit03<(5+0.006*pho_et)&&pho_IsoHcalRechit03<(3+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(3+0.001*pho_et)&&pho_IsoEcalRechit03<(5+0.006*pho_et)&&pho_IsoHcalRechit03<(3+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"IsolHighEt",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"IsolHighEt",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"IsolHighEt",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "IsolHighEt", theTriggerCutIsol);


///////////////////////////////////////////////////////////////
TString theTriggerCutIsol = "pho_IsoEcalRechit03<(5.5+0.006*pho_et)&&pho_IsoHcalRechit03<(3.5+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(3.5+0.001*pho_et)&&pho_IsoEcalRechit03<(5.5+0.006*pho_et)&&pho_IsoHcalRechit03<(3.5+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"LooseIsolHighEt",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"LooseIsolHighEt",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"LooseIsolHighEt",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "LooseIsolHighEt", theTriggerCutIsol);

///////////////////////////////////////////////////////////////////////

TString theTriggerCutIsol = "pho_IsoEcalRechit03<(6+0.006*pho_et)&&pho_IsoHcalRechit03<(4+0.0025*pho_et)";
TString theTotalIdTriggerCutIsol = "pho_hasPixelSeed==0&&pho_IsoEcalRechit<4.2&&pho_IsoHcalRechit<2.2&&((pho_isEB==1&&pho_sigmaIetaIeta<0.01)||(pho_isEE==1&&pho_sigmaIetaIeta<0.03))&&pho_IsoHollowTrkCone03<(4+0.001*pho_et)&&pho_IsoEcalRechit03<(6+0.006*pho_et)&&pho_IsoHcalRechit03<(4+0.0025*pho_et)";


do_theHistos("ratio","pho_IsoHollowTrkCone/pho_IsoHollowTrkCone03",0,5,100, theBaseCut,1,"LooserIsolHighEt",theTriggerCutIsol);
do_theHistos("relaDiff","(pho_IsoHollowTrkCone-pho_IsoHollowTrkCone03)/pho_IsoHollowTrkCone",-2,2,100, theBaseCut,1,"LooserIsolHighEt",theTriggerCutIsol);

do_theHistos("ecalIso","pho_IsoHollowTrkCone",-1,20,42, theBaseCut,0,"LooserIsolHighEt",theTotalIdTriggerCutIsol);

do_AllscatPlot("scat", "pho_IsoHollowTrkCone", "pho_IsoHollowTrkCone03", 0, 15, 100,  theBaseCut, 1, "LooserIsolHighEt", theTriggerCutIsol);




}
