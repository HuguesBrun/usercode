#include "functions.h"

double DeltaR(double phi1, double phi2, double eta1, double eta2){

  double dphi=phi2-phi1;
  if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
  if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;
  double dR=sqrt(dphi*dphi+(eta2-eta1)*(eta2-eta1));

  return dR;
}

double localtrackIsolation(TRootVertex theVertice, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot,float dRoutCone ){

        float dZcut = 0.2; //LIP
        float dxyCut   = 0.1; //TIP
        float dRinCone = 0.04;
        float etaStrip = 0.015;
        float theIsoEnergy = 0;
        for (int iTrack = 0 ;  iTrack < theTracks->GetEntriesFast() ; iTrack++){
                TRootTrack *localTrack  = (TRootTrack*) theTracks->At(iTrack);
                float dZ = fabs((localTrack->vz()-theVertice.z()) - (((localTrack->vx()-theVertice.x())*localTrack->Px() + (localTrack->vy()-theVertice.y())*localTrack->Py())) / localTrack->Pt()*localTrack->Pz()/localTrack->Pt());
                if (dZ > dZcut ) continue;
                double thePtTrack = localTrack->Pt();
                if (thePtTrack < 0 ) continue;
                //float dxy = sqrt((localTrack->vx()-theBeamSpot.x())*(localTrack->vx()-theBeamSpot.x())+(localTrack->vy()-theBeamSpot.y())*(localTrack->vy()-theBeamSpot.y()));
                //float dxy = sqrt((localTrack->vx()-theVertice.x())*(localTrack->vx()-theVertice.x())+(localTrack->vy()-theVertice.y())*(localTrack->vy()-theVertice.y()));
                //float dxy = fabs(-(localTrack->vx()-theVertice.x())*localTrack->Py()+(localTrack->vy()-theVertice.y())*localTrack->Px())/localTrack->Pt();
                float dxy = fabs(-(localTrack->vx()-theBeamSpot.x())*localTrack->Py()+(localTrack->vy()-theBeamSpot.y())*localTrack->Px())/localTrack->Pt();
                if (dxy > dxyCut) continue;
                double dR = DeltaR(localTrack->Vect().Phi(), thePhoton.Phi(), localTrack->Vect().Eta(), thePhoton.Eta());
                double dEta = fabs(localTrack->Vect().Eta()-thePhoton.Eta());
                if ((dR < dRoutCone)&&(dR > dRinCone)&&(dEta>etaStrip)) {
                        theIsoEnergy += thePtTrack;
                }
        }   
        return theIsoEnergy;
}


double  calcWorstTrackIsolation(TClonesArray* theVertices, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot, TRootVertex *theWorst){
        double theWorstIso = 0;
        int theWorstIte = 0;
        double theLocalIso;    
        TRootPhoton theLocalPhoton = thePhoton;
        for (int iVertex = 0 ; iVertex < theVertices->GetEntriesFast() ; iVertex++){
                TRootVertex *theLocalVertex = (TRootVertex*) theVertices->At(iVertex);
                TVector3 theVertex(theLocalVertex->x(),theLocalVertex->y(),theLocalVertex->z());
                theLocalPhoton.setVertex(theVertex);
                theLocalIso = localtrackIsolation(*theLocalVertex, theTracks, theLocalPhoton, theBeamSpot, 0.4);
                //cout << "local iso " << theLocalIso << endl;
                if (theLocalIso > theWorstIso) {
                        theWorstIso = theLocalIso;
                        theWorstIte = iVertex;
                }   
        }   
        TRootVertex *theLocalWorst = (TRootVertex*) theVertices->At(theWorstIte);

        *theWorst = *theLocalWorst;
        return theWorstIso;
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
void doGenInfo(TRootPhoton* myphoton, TClonesArray* mcParticles, Int_t* pho_GenId, Int_t* pho_MotherId, Int_t* pho_isGenElectron, Int_t* pho_isPromptGenPho, Int_t* pho_isFromQuarkGen, Int_t* pho_isPi0Gen, Int_t* pho_isEtaGen, Int_t* pho_isRhoGen, Int_t* pho_isOmegaGen, Float_t* pho_PromptGenIsoEnergyStatus1, Float_t* pho_PromptGenIsoEnergyStatus2, float* pho_trueE, float* pho_truePx, float* pho_truePy, float* pho_truePz, float* pho_trueEta, float* pho_truePhi,double dRcone){

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
    *pho_trueE = mygenparticle->E();
    *pho_truePx = mygenparticle->Px();
    *pho_truePy = mygenparticle->Py();
    *pho_truePz = mygenparticle->Pz();
    *pho_truePhi = mygenparticle->Phi();
    *pho_trueEta = mygenparticle->Eta();
   
    *pho_GenId = mygenparticle->type();
    *pho_MotherId = mygenparticle->motherType();
    if (abs(mygenparticle->type())==11) *pho_isGenElectron = 1;
    else *pho_isGenElectron = 0;

    if (mygenparticle->type()==22 && mygenparticle->motherType()==22) *pho_isPromptGenPho = 1;
    else *pho_isPromptGenPho = 0;

    *pho_isFromQuarkGen = 0;
    if (mygenparticle->type()==22 && mygenparticle->motherType()!=22) {

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

void matchWithAnElectron(TRootPhoton *myPhoton, TClonesArray *electrons, int *isAlsoaRecoElectron, float *pho_fBrem, float *pho_momentumCorrected, float *pho_d0, float *pho_tightEleId, float *pho_eleTrkIso, float *pho_eleEcalIso, float *pho_eleHcalIso, float *pho_eleDeltaPhiIn, float *pho_eleDeltaEtaIn, float *pho_eleHoE, float *pho_eleSigmaIetaIeta, int *pho_eleMissHits, float *pho_eleDistConvPartner, float *pho_eleDcotConvPartner){
	double dRcj = 0.1;
	unsigned int isl=1000;
	int theSCphoton = myPhoton->scIndex();
	float theSCphotonEnergy = myPhoton->superCluster()->rawEnergy();
	for (unsigned int i=0 ; i < electrons->GetEntriesFast() ; i++){
		TRootElectron *theElectron = (TRootElectron*) electrons->At(i);
	if (theElectron->Pt() < 10 ) continue;
		double DR = DeltaR(myPhoton->Phi(), theElectron->Phi(), myPhoton->Eta(), theElectron->Eta());
		if (DR < dRcj) {
//			float theSCelectronEnergy = theElectron->superCluster()->rawEnergy();
//			std::cout << "photon = " << theSCphotonEnergy << " electron = " << theSCelectronEnergy << std::endl;
//			if ( fabs(theSCphotonEnergy-theSCelectronEnergy) < 0.0001 ) {
				isl = i;
//			}
		}
//		if (theElectron->scIndex() == theSCphoton) isl = i;
	}
	if (isl != 1000){
		TRootElectron *candidate = (TRootElectron*) electrons->At(isl);
		*isAlsoaRecoElectron = 1;
		*pho_fBrem = candidate->fbrem();
		*pho_momentumCorrected = 0;//candidate->momentumCorrected();
		*pho_d0 = candidate->d0();
		*pho_tightEleId = candidate->idCutBasedFixedThresholdTight(); 
		*pho_eleTrkIso = candidate->trackIso();
		*pho_eleEcalIso = candidate->ecalIso();
		*pho_eleHcalIso = candidate->hcalIso();
		*pho_eleDeltaPhiIn = candidate->deltaEtaIn();
		*pho_eleDeltaEtaIn = candidate->deltaPhiIn();
		*pho_eleHoE = candidate->hadOverEm();
		*pho_eleSigmaIetaIeta = candidate->sigmaIetaIeta();
		*pho_eleMissHits = candidate->trackMissedInnerLayers();
		*pho_eleDistConvPartner = candidate->distConvPartner();
		*pho_eleDcotConvPartner = candidate->dcotConvPartner();

	}
}

int findMatchingWithAnHLTObjet(TRootPhoton *myPhoton, TClonesArray *HLTobject, TString filterName){
	int NbHLTSize = HLTobject->GetEntriesFast();
	if (NbHLTSize==0) return 0;
//	cout << "Size = " << NbHLTSize << endl;
	float dR;
	int isAgood = 0;
	for (int i = 0 ; i < NbHLTSize ; i++){
		TRootHLTObject *theHLT = (TRootHLTObject*) HLTobject->At(i);
		dR = DeltaR(theHLT->Phi(),myPhoton->Phi(),theHLT->Eta(),myPhoton->Eta());
		if (dR < 0.3) {
			if (filterName==theHLT->hltFilter()) isAgood=1;
		}
	}
	if (isAgood == 1) return 1;
	else return 0;
}

void findConversionMCtruth(TRootPhoton *myPhoton, TClonesArray *theMCphotons, int &pho_MCisConverted, float &pho_MCconvEoverP, float &pho_MCconvMass, float &pho_MCconvCotanTheta, float &pho_MCconvVertexX, float &pho_MCconvVertexY, float &pho_MCconvVertexZ){
float dr = 0;
int theIteMin = -1000;
float theDiff;
float theMinDiff = 100000;
for (unsigned int i =0 ; i < theMCphotons->GetEntriesFast() ; i++){
	TRootMCPhoton *theMCphoton = (TRootMCPhoton*) theMCphotons->At(i);
	dr = DeltaR(myPhoton->Phi(),theMCphoton->Phi(),myPhoton->Eta(),theMCphoton->Eta());
	if (dr < 0.3){
		theDiff = fabs(theMCphoton->Pt()-myPhoton->Pt());
		if (theDiff < theMinDiff){
			theMinDiff = theDiff;
			theIteMin = i;	
		}
	}
}
if (theIteMin == -1000){
	pho_MCisConverted = 0;
	pho_MCconvEoverP = -10000;
	pho_MCconvMass = -10000;
	pho_MCconvCotanTheta = -10000;
	pho_MCconvVertexX = -10000;
	pho_MCconvVertexY = -10000;
	pho_MCconvVertexZ = -10000;
}
else {
	TRootMCPhoton *theMCphoton = (TRootMCPhoton*) theMCphotons->At(theIteMin);
	pho_MCisConverted = 1;
	pho_MCconvEoverP = theMCphoton->convEoverP();
	pho_MCconvMass = theMCphoton->convMass();
	pho_MCconvCotanTheta = theMCphoton->convDeltaCotanTheta();
	pho_MCconvVertexX = theMCphoton->conv_vx();
	pho_MCconvVertexY = theMCphoton->conv_vy();
	pho_MCconvVertexZ = theMCphoton->conv_vz();

}
}
void findTheMCelectron(TRootPhoton *myPhoton, TClonesArray *theMCelectron, float &pho_eleMCtruthBrem, int &pho_eleMCtruthNBrem){
	int theIteMin = -1000;
	float theMinDiff  = 10000;
	float theDiff;
	float dr;
	for (unsigned int i = 0 ; i < theMCelectron->GetEntriesFast() ; i++){
		TRootMCElectron *theElectron  = (TRootMCElectron*) theMCelectron->At(i);
		dr = DeltaR(myPhoton->Phi(),theElectron->Phi(),myPhoton->Eta(),theElectron->Eta());
		if (dr < 0.3){
			theDiff = fabs(theElectron->Pt()-myPhoton->Pt());
			if (theDiff < theMinDiff){
				theMinDiff = theDiff;
				theIteMin = i;	
			}
		}
	}
	if (theIteMin > 0 ) {
		TRootMCElectron *theElectron  = (TRootMCElectron*) theMCelectron->At(theIteMin);
//		cout << "photon phi = " << myPhoton->Phi() << " electron phi " << theElectron->Phi() << endl;
		std::vector<float> theBrem = theElectron->energyLoss(); 
//$		cout << "taille = " << theBrem.size() << endl;
		pho_eleMCtruthBrem = 0;
		pho_eleMCtruthNBrem = theBrem.size();
		for (int j = 0 ; j < theBrem.size() ; j++){
	//		cout << "brem loss = " << theBrem[j] << endl;
			pho_eleMCtruthBrem += theBrem[j];
		}	
	}
}

bool isMatchingWithAMuon(TRootPhoton *thePhoton, TClonesArray *theMuons){
bool caMatch = false;
	for (int i = 0 ; i < theMuons->GetEntriesFast() ; i++){
		TRootMuon *theMuon = (TRootMuon*) theMuons->At(i);
		float dR = DeltaR(theMuon->Phi(), thePhoton->Phi(), theMuon->Eta(), thePhoton->Eta());
		if (dR < 0.4) caMatch = true;
	}

return caMatch;
}

double localtrackIsolationCIC(TRootVertex theVertice, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot,float dRoutCone ){
	
	float dZcut = 1.0; //LIP
	float dxyCut   = 0.1; //TIP
	float dRinCone = 0.02;
	float etaStrip = 0.0;
	float theIsoEnergy = 0;
	for (int iTrack = 0 ;  iTrack < theTracks->GetEntriesFast() ; iTrack++){
		TRootTrack *localTrack  = (TRootTrack*) theTracks->At(iTrack);
		float dZ = fabs((localTrack->vz()-theVertice.z()) - (((localTrack->vx()-theVertice.x())*localTrack->Px() + (localTrack->vy()-theVertice.y())*localTrack->Py())) / localTrack->Pt()*localTrack->Pz()/localTrack->Pt());
		if (dZ > dZcut ) continue;
		double thePtTrack = localTrack->Pt();
		if (thePtTrack < 0 ) continue;
		//float dxy = sqrt((localTrack->vx()-theBeamSpot.x())*(localTrack->vx()-theBeamSpot.x())+(localTrack->vy()-theBeamSpot.y())*(localTrack->vy()-theBeamSpot.y()));
		//float dxy = sqrt((localTrack->vx()-theVertice.x())*(localTrack->vx()-theVertice.x())+(localTrack->vy()-theVertice.y())*(localTrack->vy()-theVertice.y()));
		float dxy = fabs(-(localTrack->vx()-theVertice.x())*localTrack->Py()+(localTrack->vy()-theVertice.y())*localTrack->Px())/localTrack->Pt();
		if (dxy > dxyCut) continue;
	  	double dR = DeltaR(localTrack->Vect().Phi(), thePhoton.Phi(), localTrack->Vect().Eta(), thePhoton.Eta());
		double dEta = fabs(localTrack->Vect().Eta()-thePhoton.Eta());
	 	if ((dR < dRoutCone)&&(dR > dRinCone)&&(dEta>etaStrip)) {
			theIsoEnergy += thePtTrack;
		}
	}
	return theIsoEnergy;
}


double  calcWorstTrackIsolationCIC(TClonesArray* theVertices, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot, TRootVertex *theWorst){
	double theWorstIso = 0;
	int theWorstIte = 0;
	double theLocalIso;	
	TRootPhoton theLocalPhoton = thePhoton;
	for (int iVertex = 0 ; iVertex < theVertices->GetEntriesFast() ; iVertex++){
		TRootVertex *theLocalVertex = (TRootVertex*) theVertices->At(iVertex);
		TVector3 theVertex(theLocalVertex->x(),theLocalVertex->y(),theLocalVertex->z());
		theLocalPhoton.setVertex(theVertex);
		theLocalIso = localtrackIsolationCIC(*theLocalVertex, theTracks, theLocalPhoton, theBeamSpot, 0.4);
		//cout << "local iso " << theLocalIso << endl;
		if (theLocalIso > theWorstIso) {
			theWorstIso = theLocalIso;
			theWorstIte = iVertex;
		}
	} 
	TRootVertex *theLocalWorst = (TRootVertex*) theVertices->At(theWorstIte);

	*theWorst = *theLocalWorst;
	return theWorstIso;
}
int findThePhoCat(TRootPhoton photon1){
int cat = -1; 
        if (photon1.isEBPho()==1){
                if (photon1.r9()>0.94) cat = 0; 
                else cat = 1;
        }   
        else {
                if (photon1.r9()>0.94) cat = 2;    
                else cat = 3;
        }   
return cat;
}
double dRtoTrack(TRootPhoton thePhoton, TClonesArray* electrons){
	float theMinDr = 99;
	float theDr;
	for (int iTrack = 0 ; iTrack < electrons->GetEntriesFast() ; iTrack++){
		TRootElectron *theLocalElectron = (TRootElectron*) electrons->At(iTrack);
		if (theLocalElectron->trackMissedInnerLayers() > 0 ) continue;
		//if (fabs(thePhoton.superCluster()->Phi() - theLocalElectron->theSCphi())>0.01) continue;
		//if (fabs(thePhoton.superCluster()->Eta() - theLocalElectron->theSCeta())>0.01) continue;
		float deltaCalo = sqrt((thePhoton.caloPosition().X()-theLocalElectron->caloPosition().X())*(thePhoton.caloPosition().X()-theLocalElectron->caloPosition().X())+(thePhoton.caloPosition().Y()-theLocalElectron->caloPosition().Y())*(thePhoton.caloPosition().Y()-theLocalElectron->caloPosition().Y())+(thePhoton.caloPosition().Z()-theLocalElectron->caloPosition().Z())*(thePhoton.caloPosition().Z()-theLocalElectron->caloPosition().Z()));
		if (deltaCalo > 0.0001 ) continue; //no matching between the electrons and the photon :( 

		theDr = sqrt(theLocalElectron->deltaEtaIn()*theLocalElectron->deltaEtaIn()+theLocalElectron->deltaPhiIn()*theLocalElectron->deltaPhiIn());
		if (theDr < theMinDr) theMinDr = theDr;
	}
//	cout << "theDr = " << theMinDr << endl;
	return theMinDr;
}
bool photonIsPassingCIC(TRootPhoton thePhoton, TClonesArray* theVertices, TClonesArray* theTracks, TRootBeamSpot theBeamSpot, TClonesArray* electrons, int *theCutStop, int cat){
float vcicST[7][4] = { 
  {3.8,     2.2,     1.77,   1.29},    // rel combIso (good vtx)                                
  {11.7,    3.4,     3.9,    1.84},    // rel combIso (bad vtx)                                                   
  {3.5,     2.2,     2.3,    1.45},    // trkIso (good vtx)                                                           
  {0.0105,  0.0097,  0.028,  0.027},  // sigma_ie_ie                                                                  
  {0.082,   0.062,   0.065,  0.048},  // H/E       
  {0.94,    0.36,    0.94,   0.32},   // R9
  {1.,      0.062,   0.97,   0.97},    // dR to trk
};

	bool isCIC = true;
	(*theCutStop) = 0;
	double varToCut[7]; // now calc the CiC variables
	TRootVertex* theBestVertex= (TRootVertex*) theVertices->At(0); // take the default vertex 
	varToCut[0] = (thePhoton.dR03IsolationEcalRecHit() + thePhoton.dR04IsolationHcalRecHit() + localtrackIsolationCIC(*theBestVertex, theTracks, thePhoton, theBeamSpot, 0.3))*50.0/thePhoton.Et(); // comb iso respect to the selected vertex
	TRootVertex theWorstVertex;
	TRootPhoton thePhotonWithWorstVertex = thePhoton;
	varToCut[1] = (thePhoton.dR04IsolationEcalRecHit() + thePhoton.dR04IsolationHcalRecHit() + calcWorstTrackIsolationCIC(theVertices, theTracks, thePhoton, theBeamSpot, &theWorstVertex))*50; // worst comb iso 
	TVector3 theWorstCoords(theWorstVertex.x(), theWorstVertex.y(), theWorstVertex.z());
	thePhotonWithWorstVertex.setVertex(theWorstCoords);
	varToCut[1] = varToCut[1]/thePhotonWithWorstVertex.Et();
	varToCut[2] = (localtrackIsolation(*theBestVertex, theTracks, thePhoton, theBeamSpot,0.3))*50.0/thePhoton.Et();// iso track calc with the selected vertex
	//varToCut[2] = varToCut[2]*50.0/thePhoton.Et();// iso track calc with the selected vertex
	varToCut[3] = thePhoton.sigmaIetaIeta(); // sigma ieta 
	varToCut[4] = thePhoton.hoe(); // HoE
	varToCut[5] = thePhoton.r9(); // R9
	varToCut[6] = dRtoTrack(thePhoton, electrons); // calc the dR

	//cout << "photon Et = " << thePhoton.Et();
/*	for (int i = 0 ; i < 7 ; i++){
		cout << " var" << i << "= " << varToCut[i];
	}
	cout << endl;*/
	//cout << " combIso(default vtx)*50/Et=" << varToCut[0] << " 		worst comb Iso *50 /Etworst vtx=" << varToCut[1] << " iso tracker(default vtx)*50/Et=" <<  varToCut[2] << "  		sig Ieta Ieta=" << varToCut[3] << " 		HoE=" << varToCut[4] << " 		R9=" << varToCut[5] << " 		dR to Track=" << varToCut[6] << endl;
	//cout << "isoHcal = " << thePhoton.dR04IsolationHcalRecHit()  << " isoEcal = " << thePhoton.dR03IsolationEcalRecHit() << " isoTracker =  " << localtrackIsolation(*theBestVertex, theTracks, thePhoton, theBeamSpot, 0.3) << "the worst iso = " << calcWorstTrackIsolation(theVertices, theTracks, thePhoton, theBeamSpot, &theWorstVertex) << "Et worst =" << thePhotonWithWorstVertex.Et() << endl;

 	for (int i = 0 ; i < 7 ; i++){  // test if pass the CiC cuts
		//cout << "i " << varToCut[i] << " " << vcicST[i][cat-1] << endl;
		if (i < 5) {if (varToCut[i] > vcicST[i][cat]) isCIC = false; else (*theCutStop)++;}
		else {if (varToCut[i] < vcicST[i][cat]) isCIC = false; else (*theCutStop)++;} 
	}
	return isCIC;
}

