#ifndef INC_FLOCALCORR
#define INC_FLOCALCORR

#include "TVector2.h"
#include "TMath.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "CondFormats/EcalObjects/interface/EcalClusterLocalContCorrParameters.h"


namespace fLocalCorr{


  inline int leftOrRightEta( DetId theId, const EcalRecHitCollection *recHits){
	EBDetId *theEBId = new EBDetId(theId);
	EBDetId *theEBIdGau= new EBDetId(theEBId->offsetBy(-1,0));
	EBDetId *theEBIdDro = new EBDetId(theEBId->offsetBy(1,0));
	DetId *theIdGau = (DetId*) theEBIdGau; 
	DetId *theIdDro = (DetId*) theEBIdDro; 
	float EGau = EcalClusterTools::recHitEnergy((*theIdGau), recHits);
	float EDro = EcalClusterTools::recHitEnergy((*theIdDro), recHits);
	if (EGau>EDro) return 1;
	else return 0;
  }
	
  inline int leftOrRightPhi( DetId theId, const EcalRecHitCollection *recHits){
        EBDetId *theEBId = new EBDetId(theId);
        EBDetId *theEBIdInt = 0;
        EBDetId *theEBIdExt = 0;
        theEBIdInt= new EBDetId(theEBId->offsetBy(0,-1));
        theEBIdExt= new EBDetId(theEBId->offsetBy(0,1));

        DetId *theIdInt = (DetId*) theEBIdInt;
        DetId *theIdExt = (DetId*) theEBIdExt;
        float Eint = EcalClusterTools::recHitEnergy((*theIdInt), recHits);
        float Eext = EcalClusterTools::recHitEnergy((*theIdExt), recHits);
        if ( Eint > Eext ) return 1;
        else return 0;
  }
  inline float logRatioEta( DetId theId, const EcalRecHitCollection *recHits, int gaucheDroite, const CaloTopology* topology){
	float leLogRatio=-10;
	float leRatio;
	float Ein, Eout, EDro, EGau;
	EBDetId *theEBId = new EBDetId(theId);
	reco::BasicCluster ZeroCluster;
	if (!(gaucheDroite)){
		EGau = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, 0, 2 , -2, 2);
		EDro = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, -2, -1 , -2, 2);
	}
	else {
		EGau = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId,  1, 2, -2, 2);
		EDro = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, -2, 0, -2, 2);
	}
	if (theEBId->ieta() > 0 ) { Ein = EGau; Eout = EDro; }
	else { Eout = EGau; Ein = EDro; }
	if (Ein != 0) leRatio = 1.0*Eout/Ein;
		else leLogRatio =-20;
	if (leRatio>0) leLogRatio = log(leRatio);
		else leLogRatio = -20;
	return leLogRatio; 
  }
  inline float logRatioPhi( DetId theId, const EcalRecHitCollection *recHits, int innerOuter, const CaloTopology* topology){
        float leLogRatio=-10;
        float leRatio;
        float Ein, Eout;
        reco::BasicCluster ZeroCluster;
        if (innerOuter){
                Ein = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, -2, 2 , -2, -1);
                Eout = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, -2, 2 , 0, 2);
        }
        else {
                Ein = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, -2, 2, -2, 0);
                Eout = EcalClusterTools::matrixEnergy(ZeroCluster, recHits, topology, theId, -2, 2 , 1, 2);
        }
//	std::cout << "Ein = " << Ein << " Eout = " << Eout << std::endl;
        if (Ein != 0) leRatio = 1.0*Eout/Ein;
                else leLogRatio =-20;
        if (leRatio>0) leLogRatio = log(leRatio);
                else leLogRatio = -20;
        return leLogRatio;
  }
  
  inline float theTourneurLocalCorrection( const reco::SuperCluster & superCluster, const CaloSubdetectorGeometry* geom, const EcalClusterLocalContCorrParameters* params, int etaPhi ){ //const edm::EventSetup& iSetup){

  //correction factor to be returned, and to be calculated in this present function:
  double correction_factor=1.;
  double fetacor=1.; //eta dependent part of the correction factor
  double fphicor=1.; //phi dependent part of the correction factor

  //********************************************************************************************************************//
  //These local containment corrections correct a photon energy for leakage outside a 5x5 crystal cluster. They  depend on the local position in the hit crystal. The local position coordinates, called later EtaCry and PhiCry in the code, are comprised between -0.5 and 0.5 and correspond to the distance between the photon supercluster position and the center of the hit crystal, expressed in number of  crystal widthes. The correction parameters (that should be filled in CalibCalorimetry/EcalTrivialCondModules/python/EcalTrivialCondRetriever_cfi.py) were calculated using simulaion and thus take into account the effect of the magnetic field. They  only apply to unconverted t CaloGeometry*hotons in the barrel, but a use for non brem electrons could be considered (not tested yet). For more details, cf the CMS internal note 2009-013 by S. Tourneur and C. Seez

  //Beware: The user should make sure it only uses this correction factor for unconverted photons (or not breming electrons)


  const reco::CaloClusterPtr & seedbclus =  superCluster.seed();
  
  //If not barrel, return 1:
  if (TMath::Abs(seedbclus->eta()) >1.4442 ) return 1.;
 
  const math::XYZPoint position_ = seedbclus->position(); 
  double Theta = -position_.theta()+0.5*TMath::Pi();
  double Eta = position_.eta();
  double Phi = TVector2::Phi_mpi_pi(position_.phi());
 
  //Calculate expected depth of the maximum shower from energy (like in PositionCalc::Calculate_Location()):
  // The parameters X0 and T0 are hardcoded here because these values were used to calculate the corrections:
  const float X0 = 0.89; const float T0 = 7.4;
  double depth = X0 * (T0 + log(seedbclus->energy()));
  std::vector< std::pair<DetId, float> > crystals_vector = seedbclus->hitsAndFractions();
  float dphimin=999.;
  float detamin=999.;
  int ietaclosest = 0;
  int iphiclosest = 0;
  for (unsigned int icry=0; icry!=crystals_vector.size(); ++icry) {    
    EBDetId crystal(crystals_vector[icry].first);
    const CaloCellGeometry* cell=geom->getGeometry(crystal);
    GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
    double EtaCentr = center_pos.eta();
    double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
    if (TMath::Abs(EtaCentr-Eta) < detamin) {
      detamin = TMath::Abs(EtaCentr-Eta); 
      ietaclosest = crystal.ieta();
    }
    if (TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)) < dphimin) {
      dphimin = TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)); 
      iphiclosest = crystal.iphi();
    }
  }
  EBDetId crystalseed(ietaclosest, iphiclosest);
  
  // Get center cell position from shower depth
  const CaloCellGeometry* cell=geom->getGeometry(crystalseed);
  GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
    //if the seed crystal is neighbourgh of a supermodule border, don't apply the phi dependent  containment corrections, but use the larger crack corrections instead.
  int iphimod20 = TMath::Abs(iphiclosest%20);
   if ( iphimod20 <=1 ) fphicor=1.;

   else{
      double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
      double PhiWidth = (TMath::Pi()/180.);
      double PhiCry = (TVector2::Phi_mpi_pi(Phi-PhiCentr))/PhiWidth;
      if (PhiCry>0.5) PhiCry=0.5;
      if (PhiCry<-0.5) PhiCry=-0.5;
       //Some flips to take into account ECAL barrel symmetries:
      if (ietaclosest<0) PhiCry *= -1.;
      
      //Fetching parameters of the polynomial (see  CMS IN-2009/013)
      double g[5];
      for (int k=0; k!=5; ++k) g[k] = (params->params())[k+5];
      
      fphicor=0.;
      for (int k=0; k!=5; ++k) fphicor += g[k]*std::pow(PhiCry,k);
   }
  //if the seed crystal is neighbourgh of a module border, don't apply the eta dependent  containment corrections, but use the larger crack corrections instead.
  int ietamod20 = TMath::Abs(ietaclosest%20);
  if (TMath::Abs(ietaclosest) >24 && (ietamod20==5 || ietamod20==6) ) fetacor = 1.;
  
  else
    {      
      double ThetaCentr = -center_pos.theta()+0.5*TMath::Pi();
      double ThetaWidth = (TMath::Pi()/180.)*TMath::Cos(ThetaCentr);
      double EtaCry = (Theta-ThetaCentr)/ThetaWidth;    
      if (EtaCry>0.5) EtaCry=0.5;
      if (EtaCry<-0.5) EtaCry=-0.5;
      //flip to take into account ECAL barrel symmetries:
      if (ietaclosest<0) EtaCry *= -1.;
      
      //Fetching parameters of the polynomial (see  CMS IN-2009/013)
      double f[5];
      for (int k=0; k!=5; ++k) f[k] = (params->params())[k];
     
      fetacor=0.;
      for (int k=0; k!=5; ++k) fetacor += f[k]*std::pow(EtaCry,k);
    }
  
  if (etaPhi)  correction_factor  = TMath::Sqrt((params->params())[10])/(fetacor);//*fphicor);
   else correction_factor  = TMath::Sqrt((params->params())[10])/(fphicor);
  
  //*********************************************************************************************************************//
  
  //return the correction factor. Use it to multiply the cluster energy.
  return correction_factor;
   }	

  inline float theTourneurCrackCorrection( const reco::SuperCluster & superCluster, const CaloSubdetectorGeometry* geom, const EcalClusterLocalContCorrParameters* params, int etaPhi ){
     //correction factor to be returned, and to be calculated in this present function:
  double correction_factor=1.;
  double fetacor=1.; //eta dependent part of the correction factor
  double fphicor=1.; //phi dependent part of the correction factor

  //********************************************************************************************************************//
  //These ECAL barrel module and supermodule border corrections correct a photon energy for leakage outside a 5x5 crystal cluster. They  depend on the local position in the hit crystal. The hit crystal needs to be at the border of a barrel module. The local position coordinates, called later EtaCry and PhiCry in the code, are comprised between -0.5 and 0.5 and correspond to the distance between the photon supercluster position and the center of the hit crystal, expressed in number of  crystal widthes. The correction parameters (that should be filled in CalibCalorimetry/EcalTrivialCondModules/python/EcalTrivialCondRetriever_cfi.py) were calculated using simulaion and thus take into account the effect of the magnetic field. They  only apply to unconverted photons in the barrel, but a use for non brem electrons could be considered (not tested yet). For more details, cf the CMS internal note 2009-013 by S. Tourneur and C. Seez

  //Beware: The user should make sure it only uses this correction factor for unconverted photons (or not breming electrons)

  const reco::CaloClusterPtr & seedbclus =  superCluster.seed();
  
  //If not barrel, return 1:
  if (TMath::Abs(seedbclus->eta()) >1.4442 ) return 1.;

  const math::XYZPoint position_ = seedbclus->position(); 
  double Theta = -position_.theta()+0.5*TMath::Pi();
  double Eta = position_.eta();
  double Phi = TVector2::Phi_mpi_pi(position_.phi());
  
  //Calculate expected depth of the maximum shower from energy (like in PositionCalc::Calculate_Location()):
  // The parameters X0 and T0 are hardcoded here because these values were used to calculate the corrections:
  const float X0 = 0.89; const float T0 = 7.4;
  double depth = X0 * (T0 + log(seedbclus->energy()));
  std::vector< std::pair<DetId, float> > crystals_vector = seedbclus->hitsAndFractions();
  float dphimin=999.;
  float detamin=999.;
  int ietaclosest = 0;
  int iphiclosest = 0;
  for (unsigned int icry=0; icry!=crystals_vector.size(); ++icry) {    
    EBDetId crystal(crystals_vector[icry].first);
    const CaloCellGeometry* cell=geom->getGeometry(crystal);
    GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
    double EtaCentr = center_pos.eta();
    double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
    if (TMath::Abs(EtaCentr-Eta) < detamin) {
      detamin = TMath::Abs(EtaCentr-Eta); 
      ietaclosest = crystal.ieta();
    }
    if (TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)) < dphimin) {
      dphimin = TMath::Abs(TVector2::Phi_mpi_pi(PhiCentr-Phi)); 
      iphiclosest = crystal.iphi();
    }
  }
  EBDetId crystalseed(ietaclosest, iphiclosest);


  // Get center cell position from shower depth
  const CaloCellGeometry* cell=geom->getGeometry(crystalseed);
  GlobalPoint center_pos = (dynamic_cast<const TruncatedPyramid*>(cell))->getPosition(depth);
  
  //if the seed crystal isn't neighbourgh of a supermodule border, don't apply the phi dependent crack corrections, but use the smaller phi dependent local containment correction instead.
  if (ietaclosest<0) iphiclosest = 361 - iphiclosest; //inversion of phi 3 degree tilt 
  int iphimod20 = iphiclosest%20;
   if ( iphimod20 >1 ) fphicor=1.;

   else{
      double PhiCentr = TVector2::Phi_mpi_pi(center_pos.phi());
      double PhiWidth = (TMath::Pi()/180.);
      double PhiCry = (TVector2::Phi_mpi_pi(Phi-PhiCentr))/PhiWidth;
      if (PhiCry>0.5) PhiCry=0.5;
      if (PhiCry<-0.5) PhiCry=-0.5;
       //flip to take into account ECAL barrel symmetries:
      if (ietaclosest<0) PhiCry *= -1.;

      //Fetching parameters of the polynomial (see  CMS IN-2009/013)
      double g[5];
      int offset = iphimod20==0 ? 
	10 //coefficients for one phi side of a SM
	: 15; //coefficients for the other side
      for (int k=0; k!=5; ++k) g[k] = (params->params())[k+offset];
      
      fphicor=0.;
      for (int k=0; k!=5; ++k) fphicor += g[k]*std::pow(PhiCry,k);
   }

   //if the seed crystal isn't neighbourgh of a module border, don't apply the eta dependent crack corrections, but use the smaller eta dependent local containment correction instead.
  int ietamod20 = ietaclosest%20;
  if (TMath::Abs(ietaclosest) <25 || (TMath::Abs(ietamod20)!=5 && TMath::Abs(ietamod20)!=6) ) fetacor = 1.;
  
  else
    {      
      double ThetaCentr = -center_pos.theta()+0.5*TMath::Pi();
      double ThetaWidth = (TMath::Pi()/180.)*TMath::Cos(ThetaCentr);
      double EtaCry = (Theta-ThetaCentr)/ThetaWidth;    
      if (EtaCry>0.5) EtaCry=0.5;
      if (EtaCry<-0.5) EtaCry=-0.5;
      //flip to take into account ECAL barrel symmetries:
      if (ietaclosest<0) EtaCry *= -1.;
      
      //Fetching parameters of the polynomial (see  CMS IN-2009/013)
      double f[5];
      int offset = TMath::Abs(ietamod20)==5 ? 
	0 //coefficients for eta side of an intermodule gap closer to the interaction point
	: 5; //coefficients for the other eta side
      for (int k=0; k!=5; ++k) f[k] = (params->params())[k+offset];
     
      fetacor=0.;
      for (int k=0; k!=5; ++k) fetacor += f[k]*std::pow(EtaCry,k); 
    }
   if (etaPhi)   correction_factor = 1./(fetacor); //*fphicor);
    else correction_factor = 1./(fphicor);
  //*********************************************************************************************************************//
  
  //return the correction factor. Use it to multiply the cluster energy.
  return correction_factor;


  }

  inline float SLocciLocalEtaCrackCorrections(int iphi, int ieta, float logEta, int gaucheDroite){
	// first decide what correction to apply
	int leftRight = 0;
	if (gaucheDroite==1) leftRight=1; else leftRight=-1; 
	if (ieta<0) leftRight=leftRight*(-1);
	int EtaCrackPosition[5] = {0, 25, 45, 65, 85};
	int cas = 0;
	for ( int i = 0 ; i< 5 ; i++){
		if ((abs(ieta)-EtaCrackPosition[i])==-1) cas = 1;
		else if  ((abs(ieta)-EtaCrackPosition[i])==2) cas = 2;
		else if  ((abs(ieta)-EtaCrackPosition[i])==0) cas = 3;
		else if  (((abs(ieta)-EtaCrackPosition[i])==1)&&(EtaCrackPosition[i]==0)) cas = 5;
		else if  ((abs(ieta)-EtaCrackPosition[i])==1) cas = 4;
	}
	int dist = 0;       
	switch (cas) {
	     case 1 :	
		if (leftRight==1) dist = 1; else dist = 2; break;
	     case 3 :
		if (leftRight==1) dist = 2; else {dist = 5; logEta = logEta*(-1);} break;
	     case 4 :
		if (leftRight==1) dist = 5; else  { dist = 2; logEta = logEta*(-1);}  break;
	     case 2 :
		if (leftRight==1) { dist = 2; logEta = logEta*(-1); } else dist = 1; break;
	     case 5 :
		if (leftRight==1) dist = 4; else { dist = 2; logEta = logEta*(-1); } break;
	     case 0 :
		dist = 1;
	}
	float c[7], l[2];
	switch (dist) {
	     case 1:
	      // if pretty far from crack
	      c[0] =  0.989517;    c[1] = 0.00217081; c[2] = 0.00454726; c[3] = -0.000862548; 
	      c[4] = -0.000685733; c[5] = 8.16407e-5; c[6] = 3.54848e-5;
	      l[0] = -2.4; l[1] = 2.6;
	      break;
	    case 2:
	      // if approaching crack
	      c[0] =  0.989191;    c[1] = -0.000298256; c[2] = 0.00334354; c[3] = -0.00011853;
	      c[4] = -0.000354052; c[5] =  4.10106e-6;  c[6] = 0.; 
	      l[0] = -2.3; l[1] = 2.7;
	      break;
	    case 3:
	      // not used for historic reasons
	      break;
	    case 4:
	     // for the central crack
	      c[0] = 0.633061; c[1] = 0.277235; c[2] = -0.0539153; c[3] = 0.;
	      c[4] = 0.;       c[5] = 0.;       c[6] = 0.;
	      l[0] = 0.; l[1] = 2.7;
	      break;
	    case 5:
	      // for the non-central cracks
	      c[0] = 0.713747; c[1] = 0.190291;    c[2] = 0.0158045; c[3] = -0.0331736;
	      c[4] = 0.003837; c[5] = 0.000676782; c[6] = 0.;
	      l[0] = -2.583; l[1] = 2.7;
	      break;
      }  
        double             x = logEta;
	if (x<l[0])      x = l[0];
	else if (x>l[1]) x = l[1];
	double correctionFactor = c[0]+x*(c[1]+x*(c[2]+x*(c[3]+x*(c[4]+x*(c[5]+x*c[6])))));
	return correctionFactor;	
	}

	inline float SLocciLocalPhiCrackCorrections(int iphi, int ieta, float logPhi, int gaucheDroite){
		int dist  = 1;
		int where = (iphi % 20);
		switch (where) {
			case 0:
			if (gaucheDroite) dist = 2; else dist = 6;
			break;
			case 1:
			logPhi = logPhi*(-1);
			if (gaucheDroite) dist = 6; else dist = 2;
			break;
			case 2:	
			if (gaucheDroite) {dist = 2; logPhi = logPhi*(-1);} else dist = 1;
			break;
			case 19:
			if (gaucheDroite) dist = 1; else dist = 2;
			break;
		}
	float c[7], l[2];	
        switch (dist) {
             case 1:
              // if pretty far from crack
              c[0] =  0.989517;    c[1] = 0.00217081; c[2] = 0.00454726; c[3] = -0.000862548;
              c[4] = -0.000685733; c[5] = 8.16407e-5; c[6] = 3.54848e-5;
              l[0] = -2.4; l[1] = 2.6;
              break;
            case 2:
              // if approaching crack
              c[0] =  0.989191;    c[1] = -0.000298256; c[2] = 0.00334354; c[3] = -0.00011853;
              c[4] = -0.000354052; c[5] =  4.10106e-6;  c[6] = 0.;
              l[0] = -2.3; l[1] = 2.7;
              break;
	    case 6:
	      // for the ?? cracks PHI
	      c[0] = 0.806107; c[1] = -0.138045;    c[2] = -0.0255722; c[3] = -0.;
	      c[4] = 0.; c[5] = 0.; c[6] = 0.;
	      l[0] = -2.7; l[1] = 2.7;
	      break;
	    }	
        double             x = logPhi;
        if (x<l[0])      x = l[0];
        else if (x>l[1]) x = l[1];
        double correctionFactor = c[0]+x*(c[1]+x*(c[2]+x*(c[3]+x*(c[4]+x*(c[5]+x*c[6])))));
        return correctionFactor;
	}
}


#endif // INC_FLOCALCORR
