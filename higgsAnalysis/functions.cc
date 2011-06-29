#include "functions.h"

double DeltaR(double phi1, double phi2, double eta1, double eta2){
	
	double dphi=phi2-phi1;
	if (dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
	if (dphi<-TMath::Pi()) dphi=-2*TMath::Pi()-dphi;
	double dR=sqrt(dphi*dphi+(eta2-eta1)*(eta2-eta1));
	
	return dR;
}

double min(double a, double b){
	if (a<b) return a;
	else return b;
}

double CosThetaStar(TLorentzVector p1, TLorentzVector p2){
        TLorentzVector p = p1 + p2;
        TVector3 theBoost = p.BoostVector();
        TVector3 bostDir;
        if ( theBoost.Mag() != 0 ) bostDir = theBoost.Unit(); // / theBoost.Mag());
        else return -1;
        p1.Boost(-theBoost);
        if (p1.Vect().Mag()!=0) return p1.Vect().Dot(bostDir) / p1.Vect().Mag();
        else return -1;
}



double localtrackIsolation(TRootVertex theVertice, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot){
	
	float dZcut = 0.2; //LIP
	float dxyCut   = 0.1; //TIP
	float dRoutCone = 0.3;
	float dRinCone = 0.04;
	float etaStrip = 0.015;
	float theIsoEnergy = 0;
//	cout << "the photon " << thePhoton.Px() << " " << thePhoton.Py() << " " << thePhoton.Pz() << endl; 
	for (int iTrack = 0 ;  iTrack < theTracks->GetEntriesFast() ; iTrack++){
		TRootTrack *localTrack  = (TRootTrack*) theTracks->At(iTrack);
//		
//		cout << "iTrack = " << iTrack << " vz = " << localTrack->vz() << endl; //" " << localTrack->primaryVertex()->z() <<  endl;
//		cout << theVertice.z() << endl;
		float dZ = fabs(localTrack->vz()-theVertice.z());
//		cout << dZ << endl;
		if (dZ > dZcut ) continue;
//		cout << "ca passe dz = " << dZ << endl;
//		cout << localTrack->Vect().Phi() << " " << localTrack->Vect().Eta() << endl;
		double thePtTrack = localTrack->Pt();
	//	cout << "the Pt " << thePtTrack << endl;
		if (thePtTrack < 0 ) continue;
		float dxy = sqrt((localTrack->vx()-theBeamSpot.x())*(localTrack->vx()-theBeamSpot.x())+(localTrack->vy()-theBeamSpot.y())*(localTrack->vy()-theBeamSpot.y()));
		if (dxy > dxyCut) continue;
		//cout << "ca passe " << thePtTrack << " dxy " << dxy << std::endl;
	  	double dR = DeltaR(localTrack->Vect().Phi(), thePhoton.Phi(), localTrack->Vect().Eta(), thePhoton.Eta());
		double dEta = fabs(localTrack->Vect().Eta()-thePhoton.Eta());
	 	if ((dR < dRoutCone)&&(dR > dRinCone)&&(dEta>etaStrip)) {
		//	std::cout << " dR " << dR << " deta " << dEta << std::endl;
			theIsoEnergy += thePtTrack;
		}
	}

	return theIsoEnergy;
}


double  calcWorstTrackIsolation(TClonesArray* theVertices, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot){
	double theWorstIso = 0;
	double theLocalIso;	
	TRootPhoton theLocalPhoton = thePhoton;
	for (int iVertex = 0 ; iVertex < theVertices->GetEntriesFast() ; iVertex++){
		TRootVertex *theLocalVertex = (TRootVertex*) theVertices->At(iVertex);
		TVector3 theVertex(theLocalVertex->x(),theLocalVertex->y(),theLocalVertex->z());
		theLocalPhoton.setVertex(theVertex);
		theLocalIso = localtrackIsolation(*theLocalVertex, theTracks, theLocalPhoton, theBeamSpot);
		//cout << "local iso " << theLocalIso << endl;
		if (theLocalIso > theWorstIso) theWorstIso = theLocalIso;
	} 
	return theWorstIso;
}

double dRtoTrack(TRootPhoton thePhoton, TClonesArray* electrons){
	float theMinDr = 99;
	float theDr;
	cout << "nb elec = " << electrons->GetEntriesFast() << endl;
	for (int iTrack = 0 ; iTrack < electrons->GetEntriesFast() ; iTrack++){
		TRootElectron *theLocalElectron = (TRootElectron*) electrons->At(iTrack);
		cout << "missing layers " << theLocalElectron->trackMissedInnerLayers() << endl;
		if (theLocalElectron->trackMissedInnerLayers() > 0 ) continue;
		theDr = sqrt(theLocalElectron->deltaEtaOut()*theLocalElectron->deltaEtaOut()+theLocalElectron->deltaPhiOut()*theLocalElectron->deltaPhiOut());
		cout << "theDr =  " << theDr << endl;
		if (theDr < theMinDr) theMinDr = theDr;
	}
	cout << "theDr = " << theMinDr << endl;
	return theMinDr;
}



bool photonIsPassingCIC(TRootPhoton thePhoton, TClonesArray* theVertices, TClonesArray* theTracks, TRootBeamSpot theBeamSpot, TClonesArray* electrons){
	bool isCIC = true;
	double varToCut[7]; // now calc the CiC variables
	varToCut[0] = thePhoton.dR03IsolationEcalRecHit() + thePhoton.dR03IsolationHcalRecHit() + thePhoton.dR03IsolationHollowTrkCone(); // comb iso respect to the selected vertex
	varToCut[1] = thePhoton.dR03IsolationEcalRecHit() + thePhoton.dR03IsolationHcalRecHit() + calcWorstTrackIsolation(theVertices, theTracks, thePhoton, theBeamSpot); // worst comb iso 
	varToCut[2] = thePhoton.dR03IsolationHollowTrkCone(); // iso track calc with the selected vertex
	varToCut[3] = thePhoton.sigmaIetaIeta(); // sigma ieta 
	varToCut[4] = thePhoton.hoe(); // HoE
	varToCut[5] = thePhoton.r9(); // R9
	varToCut[6] = dRtoTrack(thePhoton, electrons); // calc the dR

	int cat = 0; // now find the categorie
	if (thePhoton.isEBPho()==1){
		if (thePhoton.r9() >= 0.94) cat=1;
		else cat=2;
	}
	else {
		if (thePhoton.r9() >= 0.94) cat = 3;
		else cat = 4;
	}
	/*cout << "the cat =  " << cat;
	for (int i = 0 ; i < 7 ; i++){
		cout << " the var " << i << " = " << varToCut[i];
	}
	cout << endl;*/
 	for (int i = 0 ; i < 7 ; i++){  // test if pass the CiC cuts
		//cout << "i " << varToCut[i] << " " << vcicST[i][cat-1] << endl;
		if (i < 5) {if (varToCut[i] > vcicST[i][cat-1]) isCIC = false;}
		else {if (varToCut[i] < vcicST[i][cat-1]) isCIC = false;} 
	}
	return isCIC;
}
bool photonPassingPreselection(TRootPhoton *myphoton){
//////////////     preselection //////////////////////////////////////////////////
	double etaLimits[3] = {1.4442,1.566,2.5};
	double pTmin = 20;
	double sigIetaIeta[2] = {0.013,0.03};
	double HoE = 0.15;
//////////////////////////////////////////////////////////////////////////////////
	bool passing = true;
	float photonSCeta = fabs(myphoton->superCluster()->Eta());
	//cout << "photon eta " << photonSCeta << " Et = " << myphoton->Et() << " hoe=" << myphoton->hoe() << " sigeta " << myphoton->sigmaIetaIeta() << " where " <<myphoton->isEEPho() << endl;
	if (!((photonSCeta < etaLimits[0])||((photonSCeta>etaLimits[1])&&(photonSCeta<etaLimits[2])))) passing = false;
	//cout << "cut eta " << passing << endl;
	if (!(myphoton->Et() > pTmin)) passing = false;
	//cout << "cut Pt " << passing << endl;
	if (!((myphoton->isEBPho()==1&&myphoton->sigmaIetaIeta()<sigIetaIeta[0])||(myphoton->isEEPho()==1&&myphoton->sigmaIetaIeta()<sigIetaIeta[1]))) passing = false;
	//cout << "cut sigma ieta" << passing << endl;
	if (!(myphoton->hoe() < HoE)) passing = false;
	//cout << "cut HoE " << passing << endl;
	return passing;
}

int  findGenParticle(TRootPhoton *myphoton, TClonesArray *mcParticles, TRootParticle *theOutParticle){
  TRootMCParticle* mygenparticle;
  int NbMCpartInCone=0;
  double bestPtdiff=500.0;
  int igpsl=-1;
  for (int igp=0; igp<mcParticles->GetEntriesFast(); igp++) {
    mygenparticle = (TRootMCParticle*) mcParticles->At(igp);
    if (DeltaR(mygenparticle->Phi(), myphoton->Phi(), mygenparticle->Eta(), myphoton->Eta())<0.1){
      if (mygenparticle->status()==1){
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
	*theOutParticle = *mygenparticle;  
	return 1;
  }
  else return 0;
}
