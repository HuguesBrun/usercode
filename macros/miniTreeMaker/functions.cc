#include "functions.h"

double DeltaR(double phi1, double phi2, double eta1, double eta2){

  double dphi=phi2-phi1;
  if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
  if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;
  double dR=sqrt(dphi*dphi+(eta2-eta1)*(eta2-eta1));

  return dR;
}


int* InitializeHLTinfo(TRootRun* runInfos, string* ListWantedHLTnames, int nPathWanted){

  cout << "Initializing HLT info"<<endl;
  unsigned int nPaths = runInfos->nHLTPaths();
  cout << "nPaths="<< nPaths <<endl;
  cout << "nPathWanted="<<nPathWanted<<endl;
  if (nPathWanted==0) return NULL;

  int* ListHLT = new int[nPathWanted];

  //runInfos->printHLTSummary();

  cout << "Boucle"<<endl;
  for (int i = 0; i < nPathWanted ; i++){
	ListHLT[i]=-1;
  }


  for (int ipath=0; ipath<nPaths; ipath++){
//cout << runInfos->hltNames(ipath)<<" num="<<ipath<<endl;
    for (int iwanted=0; iwanted<nPathWanted; iwanted++){
      if (ListWantedHLTnames[iwanted]==runInfos->hltNames(ipath)) ListHLT[iwanted]=ipath;

    }
  }

  cout<< "Wanted HLT :"<<endl;
  for (int iwanted=0; iwanted<nPathWanted; iwanted++){
    cout << ListWantedHLTnames[iwanted]<< " num="<< ListHLT[iwanted]<<endl;
  }

  return ListHLT;
}

double GetClosestJetEMFraction(TRootPhoton* myphoton, TClonesArray* jets, double JetPt){
           double EMFraction = -1;
              double dRcj=10;
              unsigned int isl=1000;
              for (unsigned int icj=0; icj<jets->GetEntriesFast(); icj++){
                TRootJet* myjet = (TRootJet*) jets->At(icj); 
                double DR = DeltaR(myphoton->Phi(), myjet->Phi(), myphoton->Eta(), myjet->Eta());
                if (DR<dRcj && myjet->Pt()>JetPt && myjet->ecalEnergyFraction()>0.01) {
                  dRcj=DR;
                  isl=icj;
                }
              }
              if (isl!=1000) {
                TRootJet* myjet = (TRootJet*) jets->At(isl);
		EMFraction = myjet->ecalEnergyFraction();
              }
  return EMFraction;
}

double GetDeltaRClosestJet(TRootPhoton* myphoton, TClonesArray* jets, double JetPt){
           double EMFraction = -1;
              double dRcj=10;
              unsigned int isl=1000;
              for (unsigned int icj=0; icj<jets->GetEntriesFast(); icj++){
                TRootJet* myjet = (TRootJet*) jets->At(icj); 
                double DR = DeltaR(myphoton->Phi(), myjet->Phi(), myphoton->Eta(), myjet->Eta());
                if (DR<dRcj && myjet->Pt()>JetPt && myjet->ecalEnergyFraction()>0.01) {
                  dRcj=DR;
                  isl=icj;
                }
              }
  return dRcj;
}

double GetPtOverJetPt(TRootPhoton* myphoton, TClonesArray* jets, double JetPt){
              double ptoverjetpt = -1;
              double dRcj=10;
              unsigned int isl=1000;
              for (unsigned int icj=0; icj<jets->GetEntriesFast(); icj++){
                TRootJet* myjet = (TRootJet*) jets->At(icj); 
                double DR = DeltaR(myphoton->Phi(), myjet->Phi(), myphoton->Eta(), myjet->Eta());
                if (DR<dRcj && myjet->Pt()>JetPt && myjet->ecalEnergyFraction()>0.01) {
                  dRcj=DR;
                  isl=icj;
                }
              }
              if (isl!=1000) {
                TRootJet* myjet = (TRootJet*) jets->At(isl);
                if (myjet->Pt()!=0 && myjet->ecalEnergyFraction()>0.01){
                  ptoverjetpt = myphoton->Pt()/myjet->Pt();
                }
              }
  return ptoverjetpt;
}

double GetTransverseMomentumToJetDirection(TRootPhoton *myPhoton, TClonesArray* jets, double JetPt){
               double dRcj=10;
               unsigned int isl=1000;
	       double Plongi = -1;
	       double Ptrans = -1;
               for (unsigned int icj=0; icj<jets->GetEntriesFast(); icj++){
                 TRootJet* myjet = (TRootJet*) jets->At(icj);
                 double DR = DeltaR(myPhoton->Phi(), myjet->Phi(), myPhoton->Eta(), myjet->Eta());
                 if (DR<dRcj && myjet->Pt()>JetPt && myjet->ecalEnergyFraction()>0.01) {
                   dRcj=DR;
                   isl=icj;
                 }
               }
		if (isl!=1000){
			TRootJet* myJet = (TRootJet*) jets->At(isl);
			if (myJet->Vect().Mag() > 0 ) Plongi = (myPhoton->Vect()).Dot(myJet->Vect())/myJet->Vect().Mag();
			if (Plongi != -1) Ptrans = sqrt(myPhoton->Vect().Mag2()-Plongi*Plongi);
		}
		return Ptrans;
}
void doGenInfo(TRootPhoton* myphoton, TClonesArray* mcParticles, Int_t* pho_GenId, Int_t* pho_MotherId, Int_t* pho_isGenElectron, Int_t* pho_isPromptGenPho, Int_t* pho_isFromQuarkGen, Int_t* pho_isPi0Gen, Int_t* pho_isEtaGen, Int_t* pho_isRhoGen, Int_t* pho_isOmegaGen, Float_t* pho_PromptGenIsoEnergyStatus1, Float_t* pho_PromptGenIsoEnergyStatus2, double dRcone){

  double etsumStatus1 = -1;
  double etsumStatus2 = -1;

  TRootMCParticle* mygenparticle;
  int NbMCpartInCone=0;
  double bestPtdiff=500.0;
  int igpsl=-1;
  for (int igp=0; igp<mcParticles->GetEntriesFast(); igp++) {
    mygenparticle = (TRootMCParticle*) mcParticles->At(igp);
    if (DeltaR(mygenparticle->Phi(), myphoton->Phi(), mygenparticle->Eta(), myphoton->Eta())<0.1){
      if (mygenparticle->status()==1){
	//HistoMCpartStatus1InConeId->Fill(mygenparticle->type());
	NbMCpartInCone++;
	if (fabs(mygenparticle->Pt()-myphoton->Pt())<bestPtdiff){
	  bestPtdiff=fabs(mygenparticle->Pt()-myphoton->Pt());
	  igpsl=igp;
	}
      }
    }
  }
  if (igpsl!=-1){
    mygenparticle = (TRootMCParticle*) mcParticles->At(igpsl);
    
    *pho_GenId = mygenparticle->type();
    *pho_MotherId = mygenparticle->motherType();

    if (abs(mygenparticle->type())==11) *pho_isGenElectron = 1;
    else *pho_isGenElectron = 0;

    if (mygenparticle->type()==22 && mygenparticle->motherType()==22) *pho_isPromptGenPho = 1;
    else *pho_isPromptGenPho = 0;

    if (mygenparticle->type()==22 && mygenparticle->motherType()!=22) {

      *pho_isFromQuarkGen = 0;
      *pho_isPi0Gen = 0;
      *pho_isEtaGen = 0;
      *pho_isRhoGen = 0;
      *pho_isOmegaGen = 0;

      if (mygenparticle->motherType()==21 || abs(mygenparticle->motherType())==1 || abs(mygenparticle->motherType())==2 || abs(mygenparticle->motherType())==3 || abs(mygenparticle->motherType())==4 || abs(mygenparticle->motherType())==5 || abs(mygenparticle->motherType())==6 ) *pho_isFromQuarkGen = 1;
      if (mygenparticle->motherType()==111) *pho_isPi0Gen = 1;
      if (mygenparticle->motherType()==221) *pho_isEtaGen = 1;
      if (mygenparticle->motherType()==113) *pho_isRhoGen = 1;
      if (mygenparticle->motherType()==223) *pho_isOmegaGen = 1;

    }

    if (*pho_isFromQuarkGen==1 || *pho_isPromptGenPho==1){
	etsumStatus1 = 0;
	etsumStatus2 = 0;
	//Isolated ?
	double dR, dR2;
	TRootMCParticle* photon = (TRootMCParticle*) mcParticles->At(igpsl);
	for (int igp=0; igp<mcParticles->GetEntriesFast(); igp++) {
	  if (igp!=igpsl){

	    TRootMCParticle* mygenpart = (TRootMCParticle*) mcParticles->At(igp);
	    if (mygenpart->status()==1){
	      if (mygenpart->type()!=22 || (fabs(mygenpart->Pt()-photon->Pt())>0.1 && mygenpart->type()==22)){
		dR = DeltaR(photon->Phi(), mygenpart->Phi(), photon->Eta(), mygenpart->Eta());
		if (dR<dRcone){
		  etsumStatus1 += mygenpart->Et();
		}
	      }
	    }

	    if (mygenpart->status()==2){
	      if (mygenpart->type()!=22 || (fabs(mygenpart->Pt()-photon->Pt())>0.1 && mygenpart->type()==22)){
		if  (abs(mygenpart->type())>6 && mygenparticle->motherType()!=21){
		  dR2 = DeltaR(photon->Phi(), mygenpart->Phi(), photon->Eta(), mygenpart->Eta());
		  if (dR2<dRcone){
		    etsumStatus2 += mygenpart->Et();
		  }
		}
	      }
	    }	    
	    

	  }
	}

      }

    
    
  }

  *pho_PromptGenIsoEnergyStatus1 = etsumStatus1;
  *pho_PromptGenIsoEnergyStatus2 = etsumStatus2;

  return;
}
