#include "setTDRStyle.C"
#include "TMath.h"


TString FloatToString(float number){
  ostringstream oss;
  oss << number;
  return oss.str();
}
float giveTheMax(float nb1, float nb2){
  if  ( nb1 > nb2 )  return nb1;
  else return nb2;
}

float intLumi = 1143;

TFile *myFileDATA = new TFile("diphoFile_DATA.root");
TFile *myFileDiPhoBox10to25 = new TFile("diphoFile_diPhoBox10to25.root");
TFile *myFileDiPhoBox25to250 = new TFile("diphoFile_diPhoBox25to250.root");
TFile *myFileDiPhoBox250toInf = new TFile("diphoFile_DiPhoBox250toInf.root");
TFile *myFileDiPhoBorn10to25 = new TFile("diphoFile_diPhoBorn10to25.root");
TFile *myFileDiPhoBorn25to250 = new TFile("diphoFile_diPhoBorn25to250.root");
TFile *myFileDiPhoBorn250toInf = new TFile("diphoFile_DiPhoBorn250toInf.root");
TFile *myFileGammaJet = new TFile("diphoFile_GJetPt20.root");
TFile *myFileQCDEnriched30 = new TFile("diphoFile_QCDPt30to40doubleEMEnriched.root");
TFile *myFileQCDEnriched40 = new TFile("diphoFile_QCDPt40doubleEMEnriched.root");
TFile *myFileWtoENu = new TFile("diphoFile_WtoENu.root");
TFile *myFiledrellYann = new TFile("diphoFile_DY.root");
TFile *myFileHiggsGG = new TFile("diphoFile_GluGluToHToGGM115.root");
TFile *myFileHiggsWHZHH = new TFile("diphoFile_WHZHHToGGM115.root");
TFile *myFileHiggsTTHH = new TFile("diphoFile_TTHHToGGM115.root");
TFile *myFileHiggsVBFH = new TFile("diphoFile_VBFHToGGM115.root");






int plotTheHisto(TString nomHisto, TString theXtitle, TString theYtitle){
	cout << "on va plotter " << nomHisto << endl;
	TString theNom = nomHisto;
	TString theNomData = theNom; //+"_1";
        TH1F *histoDATA = (TH1F*) myFileDATA->Get(theNomData);   	
	int NbDATA = histoDATA->Integral();
	float sigmaDATA = sqrt(NbDATA);
//File *myFile;
	float weights[15] = {
		358.20/528400, //diphoton Box
		12.37/518288,
		0.000208/515028,
		236.40/505456, // diphoton Born
		22.37/532864, 
		0.008072/525509,
		501.15/6757937, // gammaJet
		10868.00/6074670, // QCD double EM enreiched
		43571.00/39387002, 
		7899.00/2224253, //W->e nu
		1300.00/2224253, //Drell yann
		0.038617/109989, // Higgs GLuon fusion
		0.002482/102290, // Higgs strallung
		0.000236/22000, // Higgs TT
		0.002837/109828 // Higgs VBF
    };  
	TH1F *histoPrompt;
	TH1F *histoISRFSR;
	TH1F *histofake;
	TH1F *histoHiggsLocal;
	TH1F *histoDY;
	TFile *myFile;
	float NTwoReal = 0;
	float NOneRealOneFake = 0;
	float NTwoFake = 0;
	float NTotal = 0;
	float NHiggs = 0;
	float NDY = 0; 
	float SigmaTwoReal = 0;
	float SigmaOneRealOneFake = 0;
	float SigmaTwoFake = 0;
	float SigmaHiggs = 0;
	float SigmaDY = 0;
	float SigmaTotal = 0;
	float NTwoRealLumi = 0;
	float NOneRealOneFakeLumi = 0;
	float NTwoFakeLumi = 0;
	float NDYLumi = 0;
	float NTotalLumi = 0;
	float NHiggsLumi= 0;
	float sigmaTwoReal, sigmaOneRealOneFake, sigmaTwoFake, sigmaTotal, sigmaHiggs sigmaDY;
	for (int i = 0 ; i < 15 ; i ++){
	//cout << "i = " << i << endl;
	//	if (i==2) continue;
		switch (i) { 
	
		case 0 : 
			myFile = myFileDiPhoBox10to25;
			break;
		case 1 : 
			myFile = myFileDiPhoBox25to250;
			break;
		case 2 : 
			myFile = myFileDiPhoBox250toInf;
			break;
		case 3 : 
			myFile = myFileDiPhoBorn10to25;
			break;
		case 4 : 
			myFile = myFileDiPhoBorn25to250;
			break;
		case 5 : 
			myFile = myFileDiPhoBorn250toInf;
			break;
		case 6 : 
			myFile = myFileGammaJet;
			break;
		case 7 :
			myFile = myFileQCDEnriched30;
			break;
		case 8 : 
			myFile = myFileQCDEnriched40;
			break;
		case 9 :
			myFile = myFileWtoENu;
			break;
		case 10 : 
			myFile = myFiledrellYann;
			break;
		case 11 :
			myFile = myFileHiggsGG;
			break;
		case 12 :
			myFile = myFileHiggsWHZHH;
			break;
		case 13 :
			myFile = myFileHiggsTTHH;
			break;
		case 14 :
			myFile = myFileHiggsVBFH;
			break;

		}
		TH1F *histoHiggsLocal;  TH1F *histoISRFSRLocal; TString theNomfake;
		TH1F *histoDYLocal0; TH1F *histoDYLocal1; TH1F *histoDYLocal2;
		TH1F *histoPromptLocal;
		if (i>10) {
			cout << "coucou le nom " << theNom << "i " << i<< endl;
			TString theNomPrompt = theNom + "_TwoReal";
			histoHiggsLocal = (TH1F*) myFile->Get(theNomPrompt);
			NHiggs += histoHiggsLocal->Integral()*weights[i];
			SigmaHiggs += histoHiggsLocal->Integral()*weights[i]*weights[i];
			cout << "N Higgs " << NHiggs << endl;
		}
		else if  (i==10){
			cout << "coucou le nom " << theNom << "i " << i<< endl;
			TString theNomPrompt = theNom + "_TwoReal";
			
			histoDYLocal0 = (TH1F*) myFile->Get(theNomPrompt);
			TString theNomISRFSR = theNom + "_OnePromptOneFake";
			histoDYLocal1 = (TH1F*) myFile->Get(theNomISRFSR);
			TString theNomfake = theNom + "_TwoFake";
			histoDYLocal2 = (TH1F*) myFile->Get(theNomfake);
			NDY += (histoDYLocal0->Integral()+histoDYLocal1->Integral()+histoDYLocal2->Integral())*weights[i];
			cout << (histoDYLocal0->Integral()+histoDYLocal1->Integral()+histoDYLocal2->Integral()) << endl;
			cout << "NDY = " << NDY << endl;
			SigmaDY += (histoDYLocal0->Integral()+histoDYLocal1->Integral()+histoDYLocal2->Integral())*weights[i]*weights[i];
			cout << "apres" << endl;
		}
		else {
			TString theNomPrompt = theNom + "_TwoReal";
			cout << "coucou le nom " << theNom << "i " << i<< endl;
	        histoPromptLocal  = (TH1F*) myFile->Get(theNomPrompt);
			cout << "deux vrais " <<  histoPromptLocal->Integral()*weights[i] << endl;
		NTwoReal += histoPromptLocal->Integral()*weights[i];
		SigmaTwoReal += histoPromptLocal->Integral()*weights[i]*weights[i];
	        TString theNomISRFSR = theNom + "_OnePromptOneFake";
	    histoISRFSRLocal  = (TH1F*) myFile->Get(theNomISRFSR);
		NOneRealOneFake += histoISRFSRLocal->Integral()*weights[i];
		SigmaOneRealOneFake += histoISRFSRLocal->Integral()*weights[i]*weights[i];
	        TString theNomfake = theNom + "_TwoFake";
		histofakeLocal    = (TH1F*) myFile->Get(theNomfake);
		NTwoFake += histofakeLocal->Integral()*weights[i];
		SigmaTwoFake += histofakeLocal->Integral()*weights[i]*weights[i];
		}
		if (i==0) { 
			histoPrompt = new TH1F("histoPrompt","",histoPromptLocal->GetXaxis()->GetNbins(),histoPromptLocal->GetXaxis()->GetXmin(),histoPromptLocal->GetXaxis()->GetXmax());
			histoISRFSR = new TH1F("histoISRFSR","",histoPromptLocal->GetXaxis()->GetNbins(),histoPromptLocal->GetXaxis()->GetXmin(),histoPromptLocal->GetXaxis()->GetXmax());
			histofake = new TH1F("histofakeLocal","",histoPromptLocal->GetXaxis()->GetNbins(),histoPromptLocal->GetXaxis()->GetXmin(),histoPromptLocal->GetXaxis()->GetXmax());
			histoHiggs = new TH1F("histoHiggs","",histoPromptLocal->GetXaxis()->GetNbins(),histoPromptLocal->GetXaxis()->GetXmin(),histoPromptLocal->GetXaxis()->GetXmax());
			histoDY = new TH1F("histoDY","",histoPromptLocal->GetXaxis()->GetNbins(),histoPromptLocal->GetXaxis()->GetXmin(),histoPromptLocal->GetXaxis()->GetXmax());

		}
		
 	//	cout << "sous somme = " << (histoPromptLocal->Integral()+histoISRFSRLocal->Integral()+histofakeLocal->Integral()) << endl;
		
		if (i>10) histoHiggs->Add(histoHiggsLocal,weights[i]);
		else if (i == 10 ){ histoDY->Add(histoDYLocal0,weights[i]); histoDY->Add(histoDYLocal1,weights[i]); histoDY->Add(histoDYLocal2,weights[i]);}
		else {histoPrompt->Add(histoPromptLocal,weights[i]); histoISRFSR->Add(histoISRFSRLocal,weights[i]); histofake->Add(histofakeLocal,weights[i]);}
		delete histoPromptLocal; delete histoISRFSRLocal; delete histofakeLocal; delete histoHiggsLocal; delete histoDYLocal0; delete histoDYLocal1; delete histoDYLocal2;
//		 delete myFile;
	}	
	NTotal = NTwoReal + NOneRealOneFake + NTwoFake + NHiggs+ NDY;

	NTwoRealLumi = NTwoReal*intLumi*1.33;
	NOneRealOneFakeLumi = NOneRealOneFake*intLumi*1.33;
	NTwoFakeLumi = NTwoFake*intLumi;
	NHiggsLumi = NHiggs*intLumi;
	NDYLumi = NDY*intLumi*1.15;
	NTotalLumi = NTwoRealLumi + NOneRealOneFakeLumi + NTwoFakeLumi + NDYLumi;
	
	sigmaTwoRealLumi = sqrt(SigmaTwoReal)*intLumi;
	sigmaOneRealOneFakeLumi = sqrt(SigmaOneRealOneFake)*intLumi;
	sigmaTwoFakeLumi = sqrt(SigmaTwoFake)*intLumi;
	sigmaHiggsLumi = sqrt(SigmaHiggs)*intLumi;
	sigmaDYLumi = sqrt(SigmaDY)*intLumi;
	sigmaTotalLumi = sqrt(SigmaTwoReal+SigmaOneRealOneFake+SigmaTwoFake+SigmaDY)*intLumi;
	
	cout << "DATA = " << NbDATA << " +- " << sigmaDATA << endl;
	cout << " two real " << NTwoRealLumi << " +- " << sigmaTwoRealLumi <<endl;
	cout << "1 real 1 fake " << NOneRealOneFakeLumi << "+-" << sigmaOneRealOneFakeLumi << endl;
	cout << "2 fake " << NTwoFakeLumi << " +- " << sigmaTwoFakeLumi << endl;
	cout << "Higgs" << NHiggsLumi << "+-" << sigmaHiggsLumi << endl;
	cout << " total = " << NTotalLumi << " +- " << sigmaTotalLumi  << endl;
	
	TString dataPart1 = "N_{#gamma DATA} = " +FloatToString(NbDATA);
	TString dataPart2 = "  #pm " + FloatToString(sigmaDATA);
	TString dataLine = dataPart1 + dataPart2;
	
	TString TwoRealLumiPart1 = "N_{two real #gamma} = " + FloatToString(NTwoRealLumi);
	TString TwoRealLumiPart2 = " #pm " + FloatToString(sigmaTwoRealLumi);
	TString TwoRealLumiLine = TwoRealLumiPart1 + TwoRealLumiPart2;
	
	TString OneRealOneFakeLumiPart1 = "N_{one real one fake #gamma} = " + FloatToString(NOneRealOneFakeLumi);
	TString OneRealOneFakeLumiPart2 = " #pm " + FloatToString(sigmaOneRealOneFakeLumi);
	TString OneRealOneFakeLumiLine = OneRealOneFakeLumiPart1 + OneRealOneFakeLumiPart2;
	
	TString TwoFakeLumiPart1 = "N_{two fake #gamma} = " + FloatToString(NTwoFakeLumi);
	TString TwoFakeLumiPart2 = " #pm " + FloatToString(sigmaTwoFakeLumi);
	TString TwoFakeLumiLine = TwoFakeLumiPart1 + TwoFakeLumiPart2;
	
	TString HiggsLumiPart1 = "N_{Higgs} = " + FloatToString(NHiggsLumi);
	TString HiggsLumiPart2 = " #pm " + FloatToString(sigmaHiggsLumi);
	TString HiggsLumiLine = HiggsLumiPart1 + HiggsLumiPart2;
	
	TString DYLumiPart1 = "N_{DY} = " + FloatToString(NDYLumi);
	TString DYLumiPart2 = " #pm " + FloatToString(sigmaDYLumi);
	TString DYLumiLine = DYLumiPart1 + DYLumiPart2;
	
	TString TotalLumiPart1 = "N_{total MC } = " + FloatToString(NTotalLumi);
	TString TotalLumiPart2 = " #pm " + FloatToString(sigmaTotalLumi);
	TString TotalLumiLine =TotalLumiPart1 + TotalLumiPart2;
	
	      float coeff = histoDATA->Integral() / (histoPrompt->Integral() + histoISRFSR->Integral() + histofake->Integral() );

	/// here put the 
	    cout << "coeff : " << coeff << endl;
        histoPrompt->Scale(intLumi*1.3);
        histoISRFSR->Scale(intLumi*1.3);
        histofake->Scale(intLumi);
		histoHiggs->Scale(intLumi*1);
		histoDY->Scale(intLumi*1.15);

        TH1F *bas = new TH1F("bas","",histoPrompt->GetXaxis()->GetNbins(),histoPrompt->GetXaxis()->GetXmin(),histoPrompt->GetXaxis()->GetXmax());
        TH1F *milieu = new TH1F("milieu","",histoPrompt->GetXaxis()->GetNbins(),histoPrompt->GetXaxis()->GetXmin(),histoPrompt->GetXaxis()->GetXmax());
        TH1F *haut = new TH1F("haut","",histoPrompt->GetXaxis()->GetNbins(),histoPrompt->GetXaxis()->GetXmin(),histoPrompt->GetXaxis()->GetXmax());
		TH1F *hautMed = new TH1F("hautMed","",histoPrompt->GetXaxis()->GetNbins(),histoPrompt->GetXaxis()->GetXmin(),histoPrompt->GetXaxis()->GetXmax());
		TH1F *hautPlus = new TH1F("hautPlus","",histoPrompt->GetXaxis()->GetNbins(),histoPrompt->GetXaxis()->GetXmin(),histoPrompt->GetXaxis()->GetXmax());
	

	
	
	bas->Add(histoPrompt);
	milieu->Add(histoPrompt);
	haut->Add(histoPrompt);
	hautMed->Add(histoPrompt);
	hautPlus->Add(histoPrompt);
	

        milieu->Add(histoISRFSR);
        haut->Add(histoISRFSR);
        hautMed->Add(histoISRFSR);
		hautPlus->Add(histoISRFSR);

        haut->Add(histofake);
		hautMed->Add(histofake);
		hautPlus->Add(histofake);

	hautMed->Add(histoDY);
	hautPlus->Add(histoDY);
	
	
	hautPlus->Add(histoHiggs);
	
		TH1F *histoRatio = (TH1F*) histoDATA->Clone("ratio");
	histoRatio->Reset();
	histoRatio->Sumw2();
	histoRatio->Divide(histoDATA,hautMed, 1.0, 1.0);
	
	////////////////////////// Maintenant  on plot !!!!  /////////////////////////////////////
	TCanvas *c0 = new TCanvas("c0","coucou",600,600);

	    TPad *pad =new TPad("haut","haut",0,0.25,1,1);
	    pad->SetNumber(1);
	    cout << pad->GetBottomMargin() << endl;
	    pad->SetBottomMargin(0);
	    pad->Draw();
	    
	    TPad *pad2 =new TPad("bas","bas",0,0,1,0.25);
	    pad2->SetNumber(2);
	    pad2->SetTopMargin(0);
	   pad2->SetBottomMargin(0.3);
	    pad2->Draw();
	c0->cd(1);
        TLatex latexLabel;
        latexLabel.SetTextSize(0.04);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.26, 0.87, "CMS Preliminary 2010");
        latexLabel.DrawLatex(0.26, 0.82, "#sqrt{s} = 7 TeV");
        latexLabel.DrawLatex(0.26, 0.78, "L = 74 nb^{-1}");

	cout << histoDATA->GetMaximum() << endl;
	histoDATA->SetMinimum(0);
	histoDATA->SetMaximum(giveTheMax(histoDATA->GetMaximum(),haut->GetMaximum())*1.0*4/2);
//	histoDATA->GetXaxis()->SetTitle(theXtitle);
	histoDATA->GetYaxis()->SetTitle(theYtitle);
	histoDATA->SetTitleOffset(1.5,"Y");
//	histoDATA->SetTitleOffset(1,"X");
//	histoDATA->SetLabelSize(0.04,"X");
	histoDATA->SetLabelSize(0.03,"Y");
	histoDATA->SetMarkerSize(0.5);
	histoDATA->SetMarkerStyle(20);
	histoDATA->Draw("PE1");
		hautPlus->Draw("same");
		hautPlus->SetFillColor(kGreen-9);
		hautMed->Draw("same");
	    hautMed->SetFillColor(kYellow);
        haut->SetFillColor(628);
        haut->Draw("same");
        milieu->SetFillColor(kRed-7);
        milieu->Draw("same");
        bas->SetFillColor(kBlue-3);
        bas->Draw("same");

        TLatex latexLabel;
        latexLabel.SetTextSize(0.04);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.87, "CMS Preliminary 2011");
        latexLabel.DrawLatex(0.16, 0.82, "#sqrt{s} = 7 TeV");
        latexLabel.DrawLatex(0.16, 0.78, "L =  1.143 pb^{-1}");
/*	if (part==1) latexLabel.DrawLatex(0.16, 0.74, "Barrel");
	if (part==2) latexLabel.DrawLatex(0.16, 0.74, "Endcap");*/

        TLegend *lg = new TLegend(0.55, 0.7, 0.9, 0.905);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
        lg->SetShadowColor(kWhite);
        lg->AddEntry(histoDATA, "Data", "lp");
        lg->AddEntry(bas, " MC 2 reals #gamma", "f");
        lg->AddEntry(milieu, " MC one true one fake", "f");
        lg->AddEntry(haut, " MC two fakes", "f");
		lg->AddEntry(hautMed, "Drell-Yann","f");
		lg->AddEntry(hautPlus, "Higgs @ 115 GeV","f");
        lg->Draw();

       histoDATA->Draw("PE1:same:axis");
       histoDATA->Draw("PE1:same");
       
	TLatex *tl = new TLatex();
	tl->SetTextSize(0.03);
	tl->SetNDC(1);
	tl->DrawLatex(0.2,0.7,dataLine);
	tl->DrawLatex(0.2,0.67,TwoRealLumiLine);
	tl->DrawLatex(0.2,0.64,OneRealOneFakeLumiLine);
	tl->DrawLatex(0.2,0.61, TwoFakeLumiLine);
	tl->DrawLatex(0.2,0.58, DYLumiLine);
	tl->DrawLatex(0.2,0.55,HiggsLumiLine);
	tl->DrawLatex(0.2,0.52,TotalLumiLine);
      

  c0->cd(2);
       histoRatio->SetMinimum(0.0);
       histoRatio->SetMaximum(2);
       histoRatio->GetXaxis()->SetTitle(theXtitle);
       histoRatio->SetTitleOffset(1,"X"); 
       histoRatio->SetTitleSize(0.13,"X");
       histoRatio->SetLabelSize(0.11,"X"); 
       histoRatio->GetYaxis()->SetTitle("DATA / MC");
       histoRatio->SetTitleOffset(0.35,"Y"); 
       histoRatio->SetTitleSize(0.11,"Y");
       histoRatio->SetLabelSize(0.1,"Y"); 
       histoRatio->GetYaxis()->CenterTitle();
       histoRatio->SetNdivisions(509 ,"Y");
       histoRatio->Draw("EP");

       TLine *l = new TLine(histoRatio->GetXaxis()->GetXmin(),1.,histoRatio->GetXaxis()->GetXmax(),1.);
	l->SetLineWidth(1.5); 
	l->Draw("same");
	
	cout << "Nb monte Carlo = " << haut->Integral() << endl;
       
       c0->Print("gifLumi/"+theNom+"_lin.gif");
     //   c0->Print("epsLumi/"+theNom+"_lin.eps");
	
	

	c0->cd(1);
         histoDATA->SetMinimum(0.02);
         histoDATA->SetMaximum(histoDATA->GetMaximum()*1.0*2/2*exp(1.0*15/2));
	 
        c0->cd(1)->SetLogy();             
	histoDATA->Draw("PE1");
    hautPlus->Draw("same");
    hautMed->Draw("same");
	haut->Draw("same");
	milieu->Draw("same");
	bas->Draw("same");
	latexLabel.DrawLatex(0.16, 0.87, "CMS Preliminary 2011");
    latexLabel.DrawLatex(0.16, 0.82, "#sqrt{s} = 7 TeV");
    latexLabel.DrawLatex(0.16, 0.78, "L = 1.143 fb^{-1}");
    lg->Draw();
    histoDATA->Draw("PE1:same:axis");
    histoDATA->Draw("PE1:same");
	
	
	c0->cd(2);
	histoRatio->SetMaximum(2.2);
	histoRatio->Draw("same:axis");
        c0->Print("gifLumi/"+theNom+"_log.gif");                                  
        //c0->Print("epsLumi/"+theNom+"_log.eps");     
		
}
int do_dataMC_Lumi(){
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gROOT->ProcessLine(".x setTDRStyle.C");
	
//	plotTheHisto("dipho_mgg_looseCut0","M_{#gamma #gamma} (GeV) (no #sigma_{i#eta} cut)","nb #gamma / 15 GeV");

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

TString name[29] = {
"M_{#gamma #gamma} (GeV)",
"q_{t} (GeV)",
"cos(#theta *)",
"#eta *",
"lead #gamma #eta",
"trail #gamma #eta",
"lead #gamma p_{t} (GeV)",
"trail #gamma p_{t} (GeV)",
"lead #gamma R_{9}",
"trail #gamma R_{9}",
"lead #gamma NN shape output",
"trail #gamma NN shape output",
"min NN output",
"lead #gamma #sigma_{#eta,#phi}^{2}",
"trail #gamma #sigma_{#eta,#phi}^{2}",
"lead #eta - width",
"trail #eta - width",
"lead E_{max}/E_{3x3}",
"trail E_{max}/E_{3x3}",
"lead #phi - width / #eta - width",
"trail #phi - width / #eta - width",
"lead E_{max}/E_{SCraw}",
"trail E_{max}/E_{SCraw}",
"lead E_{2x2}/E_{5x5}",
"trail E_{2x2}/E_{5x5}",
"lead #lamdba^{-} / #lamdba^{+}",
"trail #lamdba^{-} / #lamdba^{+}",
"lead #lamdba^{-} / #sigma_{#eta,#eta}^{2}",
"trail #lamdba^{-} / #sigma_{#eta,#eta}^{2}"
}

TString yAxis[29] = {
"Events/5GeV",
"Events/5GeV",
"Events/0.05",
"Events/0.24",
"Events/0.2",
"Events/0.2",
"Events/5GeV",
"Events/5GeV",
"Events/0.03",
"Events/0.03",
"Events/0.05",
"Events/0.05",
"Events/0.05",
"Events/0.00002",
"Events/0.00002",
"Events/0.02",
"Events/0.02",
"Events/0.2",
"Events/0.2",
"Events/0.02",
"Events/0.02",
"Events/0.01",
"Events/0.01",
"Events/0.2",
"Events/0.2",
"Events/0.01",
"Events/0.01",
"Events/0.01",
"Events/0.01"
}

for (int i = 6 ; i < 7 ; i++){
	for (int j = 0 ; j < nbVars ; j++){
		cout << "on va ploter" << varsToPlot[j] << endl;
	//		plotTheHisto(varsToPlot[j]+Form("_loose1Cut%i",i),name[j],yAxis[j]);
		
	}	
}


/*	for (int i = 16 ; i < 17 ; i++){
		for (int j = 0 ; j < nbVars ; j++){
			cout << "on va ploter" << varsToPlot[j] << endl;
			plotTheHisto(varsToPlot[j]+Form("_romaRhoCut%i",i),name[j],yAxis[j]);
			
		}	
	}*/

	plotTheHisto("vertex","nb of vertex","Events");
	
	plotTheHisto("higgs1_CIC","M_{#gamma #gamma} (GeV)","Events / GeV");
	plotTheHisto("higgs2_CIC","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs3_CIC","M_{#gamma #gamma} (GeV)","Events / 3 GeV");
	plotTheHisto("higgs4_CIC","M_{#gamma #gamma} (GeV)","Events / 4 GeV");
	plotTheHisto("higgs5_CIC","M_{#gamma #gamma} (GeV)","Events / 5 GeV");
	
	plotTheHisto("higgs2_CIC_cat0","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2_CIC_cat1","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2_CIC_cat2","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2_CIC_cat3","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
    
    plotTheHisto("higgs2_CICL_cat0","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2_CICL_cat1","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2_CICL_cat2","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2_CICL_cat3","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	
	plotTheHisto("higgs2NN_CIC_cat0","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2NN_CIC_cat1","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2NN_CIC_cat2","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs2NN_CIC_cat3","M_{#gamma #gamma} (GeV)","Events / 2 GeV");
    
    plotTheHisto("higgs2NN_CIC","M_{#gamma #gamma} (GeV)","Events / 2 GeV");


    plotTheHisto("minNN_CIC_cat0","min NN output", "Events / 0.05 GeV");
    plotTheHisto("minNN_CIC_cat1","min NN output", "Events / 0.05 GeV");
    plotTheHisto("minNN_CIC_cat2","min NN output", "Events / 0.05 GeV");
    plotTheHisto("minNN_CIC_cat3","min NN output", "Events / 0.05 GeV");

	
/*	plotTheHisto("higgs1_romaRho","M_{#gamma #gamma} (GeV)","Events / GeV");
	plotTheHisto("higgs2_romaRho","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgs3_romaRho","M_{#gamma #gamma} (GeV)","Events / 3 GeV");
	plotTheHisto("higgs4_romaRho","M_{#gamma #gamma} (GeV)","Events / 4 GeV");
	plotTheHisto("higgs5_romaRho","M_{#gamma #gamma} (GeV)","Events / 5 GeV");

	plotTheHisto("higgsNN1_CIC","M_{#gamma #gamma} (GeV)","Events / GeV");
	plotTheHisto("higgsNN2_loose1","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgsNN3_loose1","M_{#gamma #gamma} (GeV)","Events / 3 GeV");
	plotTheHisto("higgsNN4_loose1","M_{#gamma #gamma} (GeV)","Events / 4 GeV");
	plotTheHisto("higgsNN5_loose1","M_{#gamma #gamma} (GeV)","Events / 5 GeV");
	
	plotTheHisto("higgsNN1_romaRho","M_{#gamma #gamma} (GeV)","Events / GeV");
	plotTheHisto("higgsNN2_romaRho","M_{#gamma #gamma} (GeV)","Events / 2 GeV");	
	plotTheHisto("higgsNN3_romaRho","M_{#gamma #gamma} (GeV)","Events / 3 GeV");
	plotTheHisto("higgsNN4_romaRho","M_{#gamma #gamma} (GeV)","Events / 4 GeV");
	plotTheHisto("higgsNN5_romaRho","M_{#gamma #gamma} (GeV)","Events / 5 GeV");*/

}