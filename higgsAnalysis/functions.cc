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



double localtrackIsolation(TRootVertex theVertice, TClonesArray* theTracks, TRootPhoton thePhoton, TRootBeamSpot theBeamSpot,float dRoutCone ){
	
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

int findTheDiphoCat(TRootPhoton photon1, TRootPhoton photon2){
int cat = -1;
	if ((photon1.isEBPho()==1)&&(photon2.isEBPho()==1)){
		if (min(photon1.r9(),photon2.r9())>0.94) cat = 0;	
		else cat = 1;
	}
	else {
		if (min(photon1.r9(),photon2.r9())>0.94) cat = 2;	
		else cat = 3;
	}
return cat;
}

bool photonIsPassingCIC(TRootPhoton thePhoton, TClonesArray* theVertices, TClonesArray* theTracks, TRootBeamSpot theBeamSpot, TClonesArray* electrons, int *theCutStop, int cat){
	bool isCIC = true;
	(*theCutStop) = 0;
	double varToCut[7]; // now calc the CiC variables
	TRootVertex* theBestVertex= (TRootVertex*) theVertices->At(0); // take the default vertex 
	varToCut[0] = (thePhoton.dR03IsolationEcalRecHit() + thePhoton.dR04IsolationHcalRecHit() + localtrackIsolation(*theBestVertex, theTracks, thePhoton, theBeamSpot, 0.3))*50.0/thePhoton.Et(); // comb iso respect to the selected vertex
	TRootVertex theWorstVertex;
	TRootPhoton thePhotonWithWorstVertex = thePhoton;
	varToCut[1] = (thePhoton.dR04IsolationEcalRecHit() + thePhoton.dR04IsolationHcalRecHit() + calcWorstTrackIsolation(theVertices, theTracks, thePhoton, theBeamSpot, &theWorstVertex))*50; // worst comb iso 
	TVector3 theWorstCoords(theWorstVertex.x(), theWorstVertex.y(), theWorstVertex.z());
	thePhotonWithWorstVertex.setVertex(theWorstCoords);
	varToCut[1] = varToCut[1]/thePhotonWithWorstVertex.Et();
	varToCut[2] = (localtrackIsolation(*theBestVertex, theTracks, thePhoton, theBeamSpot,0.3))*50.0/thePhoton.Et();// iso track calc with the selected vertex
	varToCut[2] = varToCut[2]*50.0/thePhoton.Et();// iso track calc with the selected vertex
	varToCut[3] = thePhoton.sigmaIetaIeta(); // sigma ieta 
	varToCut[4] = thePhoton.hoe(); // HoE
	varToCut[5] = thePhoton.r9(); // R9
	varToCut[6] = dRtoTrack(thePhoton, electrons); // calc the dR

/*	for (int i = 0 ; i < 7 ; i++){
		cout << " var" << i << "= " << varToCut[i];
	}
	cout << endl;**/

 	for (int i = 0 ; i < 7 ; i++){  // test if pass the CiC cuts
		//cout << "i " << varToCut[i] << " " << vcicST[i][cat-1] << endl;
		if (i < 5) {if (varToCut[i] > vcicST[i][cat]) isCIC = false; else (*theCutStop)++;}
		else {if (varToCut[i] < vcicST[i][cat]) isCIC = false; else (*theCutStop)++;} 
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
	if (!((photonSCeta < etaLimits[0])||((photonSCeta>etaLimits[1])&&(photonSCeta<etaLimits[2])))) passing = false;
	if (!(myphoton->Et() > pTmin)) passing = false;
	if (!((myphoton->isEBPho()==1&&myphoton->sigmaIetaIeta()<sigIetaIeta[0])||(myphoton->isEEPho()==1&&myphoton->sigmaIetaIeta()<sigIetaIeta[1]))) passing = false;
	if (!(myphoton->hoe() < HoE)) passing = false;
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
