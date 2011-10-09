#include "setTDRStyle.C"
#include "TMath.h"



using namespace RooFit;


TString FloatToString(float number){
        ostringstream oss;
        oss << number;
        return oss.str();
}


void do_likelihoodFun(TString VarNom, TString LatexNom, float Xmin, float Xmax, int bins, int part){

        // create the x variable (we will use it for s4ratio)
	RooRealVar x("x","x", Xmin, Xmax);


        //declare and fill the chain for Summmer11 
	TChain *chain = new TChain("photons");
	chain->Add("/sps/cms/jfan/CMSSW_4_2_3/src/TMVA/MCvsDATA/ComMC/MCDataSet/423MC/diPhoBorn25to250/output_*.root");
	/// declare and fill the chain for Spring11
	TChain *chain2 = new TChain("photons");
	chain2->Add("/sps/cms/jfan/CMSSW_4_2_3/src/TMVA/MCvsDATA/ComMC/MCDataSet/414MC/diPhoBorn25to250/output_*.root");


	//now put s4 ratio in a histo 
	TH1D *hSpring11 = new TH1D("hSpring11", "", bins, Xmin, Xmax);
        TString TheCuts = "pho_et>20 && (abs(pho_eta)<=2.5) && (!((abs(pho_eta)>1.4442)&&(abs(pho_eta)<1.566)))";
        if(part==0){
              TString Cuts = TheCuts + " && pho_isEB==1";
              TString theNom = VarNom + "_EB";
              TString theLatexNom = LatexNom + "_EB";
        }
        if(part==1){
              TString Cuts = TheCuts + " && pho_isEE==1";
              TString theNom = VarNom + "_EE";
              TString theLatexNom = LatexNom + "_EE";
        }        
        TString NomSpring11 = theNom + "_Spring11";
        TCanvas *c1 = new TCanvas(NomSpring11, "coucou", 600, 600);
        c1->cd();
        TString DrawVarNom = VarNom + ">> hSpring11";
	chain2->Draw(DrawVarNom, Cuts);
        cout << Cuts <<endl;     
        delete c1; 


        //this histo will be the pdf of the "good" shape
	RooDataHist dataHistSpring11("dataHistSpring11","the data histo spring 11",x, hSpring11);
	RooHistPdf theHistoPdf("theHistoPdf","the Histo pdf",x,dataHistSpring11);
	RooAbsReal*  nll;
	TH1D *hLL = new TH1D("h","", 40, 0.795, 1.195);
	//TH1D *hLL = new TH1D("h","", 30, 0.9895, 1.0195);
        TString NomSummer11 = theNom + "_Summer11";
        TCanvas *c2 = new TCanvas(NomSummer11, "coucou", 600, 600);
        c2->cd();
        float minNLL = 100000;
        float theGoodI;
	//for (int i = -20 ; i < 20 ; i++){	
	for (int i = -10 ; i < 20 ; i++){	
		TH1D *hSummer11 = new TH1D("hSummer11","", bins, Xmin, Xmax);
                chain->Draw(Form("%f*", 1.0+0.01*i) + VarNom + ">>hSummer11", Cuts);
                //chain->Draw(Form("%f*", 1.0+0.001*i) + VarNom + ">>hSummer11", Cuts);

		// this histo will become the "data"
		RooDataHist dataHistSummer11("dataHistSummer11","the data histo summer 11",x, hSummer11);
		//theHistoPdf.fitTo(dataHistSummer11,Extended()) ;
		nll = theHistoPdf.createNLL(dataHistSummer11);
		hLL->Fill(1.0+0.01*i, -1.0*nll->getVal());
		//hLL->Fill(1.0+0.001*i, -1.0*nll->getVal());
		delete hSummer11;
                if (minNLL > nll->getVal()) 
                {   
                        minNLL =  nll->getVal();
                        //theGoodI = 1.0+0.01*i;
                        theGoodI = 1.0+0.001*i;
                }    
                cout<< "par = " << 1.0+0.01*i << "     nll = " << -1.0*nll->getVal() <<endl;
                //cout<< "par = " << 1.0+0.001*i << "     nll = " << -1.0*nll->getVal() <<endl;



		/*TCanvas *c0 = new TCanvas("c0","coucou",600,600);
		c0->SetFillColor(0);
		RooPlot* frame1 = x.frame(Bins(40),Title("s4 ratio")) ;
		dataHistSummer11.plotOn(frame1);
		theHistoPdf.plotOn(frame1);
		frame1->Draw();*/	
	}
        delete c2;


       
        RooRealVar LLvalue("LLvalue","value of parameter", 0.895, 1.195);
        //RooRealVar LLvalue("LLvalue","value of parameter", 0.9895, 1.0195);
       /*
	RooRealVar mean("mean","mean", 0.95, 0.8, 1.1) ;
	RooRealVar w("w","width", 0.24, 0.01, 0.5 ) ;
	RooGaussian theGauss("theGauss","higgs signal",LLvalue,mean,w) ;
        */
	RooDataHist LLdataHist("LLdataHist","likelihood value",LLvalue, hLL);
	//theGauss.fitTo(LLdataHist);
        

        TCanvas *c0 = new TCanvas(theNom, theLatexNom, 600, 600);
        c0->cd();
        c0->SetFillColor(0);
        RooPlot* frame1 = LLvalue.frame(Bins(20),Title("NLL")) ;
        LLdataHist.plotOn(frame1);
        LLdataHist.plotOn(frame1, Name("LLdataHist"));
        /*
        theGauss.plotOn(frame1);
        theGauss.plotOn(frame1, Name("theGauss"));
        mean.setPlotLabel("Gaus Mean");
        theGauss.paramOn(frame1,Parameters(RooArgSet(mean)),Layout(0.3,0.65,0.55));
        float chi2 = frame1->chiSquare("theGauss", "LLdataHist", 2);
        */
        frame1->Draw();
        cout << "Good Value = " << theGoodI <<endl;

     
        TLatex latexLabel;
        latexLabel.SetTextSize(0.04);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.3, 0.6, theLatexNom);
        //TString Chi2 = "#chi^{2}/ndof = " + FloatToString(chi2);
        //latexLabel.DrawLatex(0.3, 0.4, Chi2);
        TString maxValue = "Max = " + FloatToString(theGoodI);
        latexLabel.DrawLatex(0.3, 0.5, maxValue);
      
      

        delete hSpring11;
        delete hLL;


/*                TH1D *hSummer11 = new TH1D("hSummer11","",40,0.4,1);
                chain->Draw(Form("%f*pho_e2x2/pho_e5x5>>hSummer11",1.00643),"pho_et>20 && (abs(pho_eta)<=2.5) && (!((abs(pho_eta)>1.4442)&&(abs(pho_eta)<1.566)))");
                RooDataHist dataHistSummer11("dataHistSummer11","the data histo summer 11",x, hSummer11);
                nll = theHistoPdf.createNLL(dataHistSummer11);
                delete hSummer11;
                TCanvas *c1 = new TCanvas("c1","coucou",600,600);
                c1->SetFillColor(0);
                RooPlot* frame2 = x.frame(Bins(40),Title("s4 ratio")) ;
                dataHistSummer11.plotOn(frame2);
                theHistoPdf.plotOn(frame2);
                frame2->Draw();
*/

}


void do_likelihood(){

/*
    do_likelihoodFun("pho_e2x2/pho_e5x5", "ratioS4", 0.5, 0.95, 45, 0); 
   
    do_likelihoodFun("pho_e2x2/pho_e5x5", "ratioS4", 0.5, 0.95, 45, 1 );
    do_likelihoodFun("pho_r9", "r9", 0.2, 1, 50, 0 );
    do_likelihoodFun("pho_r9", "r9", 0.2, 1, 50, 1 );
     
    //do_likelihoodFun("pho_etawidth", "etawidth", 0.005, 0.015, 20, 0 );
    //do_likelihoodFun("pho_etawidth", "etawidth", 0.007, 0.035, 28, 1 );

    do_likelihoodFun("pho_phiwidth", "phiwidth", 0.007, 0.095, 44, 0 );
    do_likelihoodFun("pho_phiwidth", "phiwidth", 0.007, 0.095, 44, 1 );
  
    //do_likelihoodFun("pho_phiwidth/pho_etawidth", "brem", 0.7, 9.7, 45, 0 );
  
    //do_likelihoodFun("pho_phiwidth/pho_etawidth", "brem", 0.7, 5.7, 50, 1 );
  */ 
    //do_likelihoodFun("pho_cPP", "cPP", 0.00005, 0.0005, 45, 0 );
    //do_likelihoodFun("pho_cPP", "cPP", 0.0001, 0.003, 29, 1 );
    //do_likelihoodFun("pho_cEE", "cEE", 0.00002, 0.0002, 36, 0 );
    //do_likelihoodFun("pho_cEE", "cEE", 0.00005, 0.0012, 10, 1 );


    //do_likelihoodFun("pho_SCbr", "SCbr", 0.7, 7.7, 35, 0 );
    //do_likelihoodFun("pho_SCbr", "SCbr", 0.6, 5.2, 23, 1 );

    //do_likelihoodFun("pho_eMax/pho_SCEraw", "ratioSeed", 0.12, 0.84, 36, 0 );
    //do_likelihoodFun("pho_eMax/pho_SCEraw", "ratioSeed", 0.12, 0.88, 37, 1 );

    do_likelihoodFun("pho_sigmaIetaIeta", "sigmaIetaIeta", 0.006, 0.012, 30, 0); 
    //do_likelihoodFun("pho_sigmaIetaIeta", "sigmaIetaIeta", 0.018, 0.039, 30, 1); 
} 























