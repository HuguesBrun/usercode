#ifndef INC_MYCALOTOOLS
#define INC_MYCALOTOOLS

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
 
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// SRP RecHits flags collection
#include "DataFormats/EcalDigi/interface/EBSrFlag.h"
#include "DataFormats/EcalDigi/interface/EESrFlag.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h" 
 
#include "CLHEP/Geometry/Transform3D.h"
 
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"


namespace myCaloTools{
  inline int deltaIphi(int iphi0, int iphi){
    int delta0 = (int) fabs(iphi0 - iphi);
    int delta1;
    if (delta0 > 180 ) {
	delta1 = 360 - delta0;
    }		
    else{
	delta1 = delta0;
    }
   return delta1;
  }
  inline int contourinfEndcap(const reco::BasicCluster &cluster, DetId id) {
        unsigned int nbint = 0;
        const EEDetId *seedXtal = new EEDetId(id.rawId());
        int iX0 = seedXtal->ix();
        int iY0 = seedXtal->iy();
  //       std::cout << iphi0 << " " << ieta0 << std::endl;
        const std::vector< std::pair<DetId, float> > listcristaux = cluster.hitsAndFractions();
           for ( unsigned int k = 0 ; k < listcristaux.size() ; k++){
                DetId leid = listcristaux[k].first;
                const EEDetId *testBC = new EEDetId(leid.rawId());
                if ((fabs(iX0-testBC->ix())<=1)&&(fabs(iY0-testBC->iy())<=1)&&(!(*testBC==*seedXtal))){
//                      std::cout << testBC->iphi() << " " << testBC->ieta() << "taux = " << listcristaux[k].second << std::endl;
                        nbint++;
                }

           }


        return nbint;
  }

  inline int contoursupEndcap(const reco::BasicCluster &cluster, DetId id) {
        unsigned int nbint = 0;
        const EEDetId *seedXtal = new EEDetId(id.rawId());
        int iX0 = seedXtal->ix();
        int iY0 = seedXtal->iy();
        //std::cout << iphi0 << " " << ieta0 << std::endl;
        const std::vector< std::pair<DetId, float> > listcristaux = cluster.hitsAndFractions();
           for ( unsigned int k = 0 ; k < listcristaux.size() ; k++){
                DetId leid = listcristaux[k].first;
                const EEDetId *testBC = new EEDetId(leid.rawId());
                if ((fabs(iX0-testBC->ix())<=2)&&(fabs(iY0-testBC->iy())<=2)&&(!((fabs(iX0-testBC->ix())<=1)&&(fabs(iY0-testBC->iy())<=1)))){
//                      std::cout << testBC->iphi() << " " << testBC->ieta() << std::endl;
                        nbint++;
                }

           }


        return nbint;

  }

  inline int contourinf(const reco::BasicCluster &cluster, DetId id,  bool isBarrel) {
	unsigned int nbint = 0;
        if (isBarrel) {
            const EBDetId *seedXtal = new EBDetId(id.rawId());
     	    int iphi0 = seedXtal->iphi();
	    int ieta0 = seedXtal->ieta();
  //        std::cout << iphi0 << " " << ieta0 << std::endl;
            const std::vector< std::pair<DetId, float> > listcristaux = cluster.hitsAndFractions();
            for ( unsigned int k = 0 ; k < listcristaux.size() ; k++){
		DetId leid = listcristaux[k].first;
                const EBDetId *testBC = new EBDetId(leid.rawId());
		if ((deltaIphi(iphi0,testBC->iphi())<=1)&&(fabs(ieta0-testBC->ieta())<=1)&&(!(*testBC==*seedXtal))){
//                      std::cout << testBC->iphi() << " " << testBC->ieta() << "taux = " << listcristaux[k].second << std::endl;
			nbint++;
		}	
 
           }
         }
	else nbint = contourinfEndcap(cluster, id);
 
	return nbint;
  }

  inline int contoursup(const reco::BasicCluster &cluster, DetId id, bool isBarrel) {
        unsigned int nbint = 0;
	if (isBarrel) {
            const EBDetId *seedXtal = new EBDetId(id.rawId());
            int iphi0 = seedXtal->iphi();
            int ieta0 = seedXtal->ieta();
         	//std::cout << iphi0 << " " << ieta0 << std::endl;
           const std::vector< std::pair<DetId, float> > listcristaux = cluster.hitsAndFractions();
           for ( unsigned int k = 0 ; k < listcristaux.size() ; k++){
                DetId leid = listcristaux[k].first;
                const EBDetId *testBC = new EBDetId(leid.rawId());
                if ((deltaIphi(iphi0,testBC->iphi())<=2)&&(fabs(ieta0-testBC->ieta())<=2)&&(!((deltaIphi(iphi0,testBC->iphi())<=1)&&(fabs(ieta0-testBC->ieta())<=1)))){
//			std::cout << testBC->iphi() << " " << testBC->ieta() << std::endl;
                        nbint++;
                }

           }
        }
	else nbint = contoursupEndcap(cluster, id);	

        return nbint;

  }

  inline float rookEnergy(const EcalRecHitCollection *recHits, const CaloTopology* topology, DetId id){
     CaloNavigator<DetId> cursor = CaloNavigator<DetId>( id, topology->getSubdetectorTopology( id ) );
     float theEnergy = 0;
     float localEnergy;
     cursor.home();
     cursor.offsetBy( 0, 1 );
     localEnergy = EcalClusterTools::recHitEnergy( *cursor, recHits );
     if (localEnergy>theEnergy) theEnergy = localEnergy;
     cursor.home();
     cursor.offsetBy( 0, -1 );
     localEnergy = EcalClusterTools::recHitEnergy( *cursor, recHits );
     if (localEnergy>theEnergy) theEnergy = localEnergy;
     cursor.home();
     cursor.offsetBy( 1, 0 );
     localEnergy = EcalClusterTools::recHitEnergy( *cursor, recHits );
     if (localEnergy>theEnergy) theEnergy = localEnergy;
     cursor.home();
     cursor.offsetBy( -1, 0 );
     localEnergy = EcalClusterTools::recHitEnergy( *cursor, recHits );
     if (localEnergy>theEnergy) theEnergy = localEnergy;
     return theEnergy; 
  }
  inline float recupChi2(const EcalRecHitCollection *recHits, DetId id){
      EcalRecHitCollection::const_iterator it;
      for (it = recHits->begin(); it != recHits->end(); it++){//tourn	
	if (id == it->id() ) break;
      }
     return it->chi2(); 
  }

  inline float recupTime(const EcalRecHitCollection *recHits, DetId id){
      EcalRecHitCollection::const_iterator it;
      for (it = recHits->begin(); it != recHits->end(); it++){//tourn	
	if (id == it->id() ) break;
      }
     return it->time(); 
  }
  inline  int recupFlag(const EcalRecHitCollection *recHits, DetId id){
      EcalRecHitCollection::const_iterator it;
      for (it = recHits->begin(); it != recHits->end(); it++){//tourn
        if (id == it->id() ) break;
      }
     return it->flags();
  }


}

#endif //INC_MYCALOTOOLS

