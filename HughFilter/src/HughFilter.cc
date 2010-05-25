// -*- C++ -*-
//
// Package:    HughFilter
// Class:      HughFilter
// 
/**\class HughFilter HughFilter.cc hugues/HughFilter/src/HughFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Hugues Louis Brun
//         Created:  Wed Dec 16 11:42:26 CET 2009
// $Id$
//
//

// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/Handle.h"



/*#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Geometry/Point3D.h"*/

#include "DataFormats/L1GlobalTrigger/interface/L1GtTechnicalTrigger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
/*#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"*/
//include "RecoEcal/EgammaClusterProducers/interface/HybridClusterProducer.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackExtra.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"


//
// class declaration
//

class HughFilter : public edm::EDFilter {
   public:
      explicit HughFilter(const edm::ParameterSet&);
      ~HughFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  
     // ----------member data ---------------------------
      edm::InputTag barrelEcalHits_;
      edm::InputTag endcapEcalHits_;
      std::string superClusterProducer_;
      std::string superClusterCollection_;
      std::string superClusterProducerEE_;
      std::string superClusterCollectionEE_;
      edm::InputTag triggerL1Tag_;
      float etmin_;
      unsigned int Xtalmin_;
      int nEvent_;	
      int runNumber_;
      bool keepOnly20_;	
      std::vector<int> eventNumbers_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HughFilter::HughFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  barrelEcalHits_ = iConfig.getParameter<edm::InputTag>("barrelEcalHits");
  endcapEcalHits_ = iConfig.getParameter<edm::InputTag>("endcapEcalHits");
  superClusterProducer_ = iConfig.getParameter<std::string>("superClusterProducer");
  superClusterCollection_ = iConfig.getParameter<std::string>("superClusterCollection");
  superClusterProducerEE_ = iConfig.getParameter<std::string>("superClusterProducerEndcap");
  superClusterCollectionEE_ = iConfig.getParameter<std::string>("superClusterCollectionEndcap");

  triggerL1Tag_ = iConfig.getParameter<edm::InputTag>("L1triggerResults");
  etmin_ = iConfig.getParameter<double>("ScetThreshold");
  Xtalmin_ = iConfig.getParameter<int>("NbXtalThreshold");
  nEvent_ = iConfig.getParameter<int>("nEvent");
  runNumber_ = iConfig.getParameter<int>("runNumber");
  eventNumbers_ = iConfig.getParameter<std::vector<int> >("eventNumbers");
  keepOnly20_ = iConfig.getParameter<bool>("keepOnlyTwenty");
}


HughFilter::~HughFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HughFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif*/
//  const reco::SuperClusterCollection superClusters; 


/*    edm::ESHandle<CaloTopology> theCaloTopology_;
    iSetup.get<CaloTopologyRecord>().get(theCaloTopology_);
  const CaloTopology *topology = theCaloTopology_.product();




   edm::Handle<EcalRecHitCollection> rhcHandleEB_;
   iEvent.getByLabel(barrelEcalHits_,rhcHandleEB_);

   if (!(rhcHandleEB_.isValid())) {
     edm::LogInfo("EnergyScaleAnalyzer") << "could not get a handle on the EcalRecHitCollection!";
   }
   const EcalRecHitCollection *EBhit_collection = rhcHandleEB_.product();



   edm::Handle<EcalRecHitCollection> rhcHandleEE_;
   iEvent.getByLabel(endcapEcalHits_,rhcHandleEE_);

  if (!(rhcHandleEE_.isValid())) {
    edm::LogInfo("EnergyScaleAnalyzer") << "could not get a handle on the EcalRecHitCollection!";
  }
  const EcalRecHitCollection *EEhit_collection = rhcHandleEE_.product();





  Handle<reco::SuperClusterCollection> pSuperClusters;
  try {
    iEvent.getByLabel(superClusterProducer_, superClusterCollection_, pSuperClusters);
  } catch ( cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer")
      << "L210 Error! can't get collection with label "
      << superClusterCollection_.c_str() ;
  }
  const reco::SuperClusterCollection* superClusters = pSuperClusters.product();


 Handle<reco::SuperClusterCollection> pSuperClustersEE;
  try {
    iEvent.getByLabel(superClusterProducerEE_, superClusterCollectionEE_, pSuperClustersEE);
  } catch ( cms::Exception& ex ) {
    edm::LogError("HughFilter")
      << "L210 Error! can't get collection with label "
      << superClusterCollectionEE_.c_str() ;
  }
  const reco::SuperClusterCollection* superClustersEE = pSuperClustersEE.product();


 

///////////////  L1 trigger //////////////////////////////////////////

    edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
    iEvent.getByLabel(triggerL1Tag_, gtReadoutRecord);
    const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
//    if (!((technicalTriggerWordBeforeMask.at(0))&&((technicalTriggerWordBeforeMask.at(40))||(technicalTriggerWordBeforeMask.at(41)))&&(!((technicalTriggerWordBeforeMask.at(36))||(technicalTriggerWordBeforeMask.at(37))||(technicalTriggerWordBeforeMask.at(38))||(technicalTriggerWordBeforeMask.at(39))))))  return false;


//     etmin_ = 10;
     bool isMoreEt = false;
     bool isGoodR9 = false;
     unsigned int nbXtals = 0;
     bool isGoodSC = false;
     
     for(reco::SuperClusterCollection::const_iterator em = superClusters->begin();
        em != superClusters->end(); ++em) {
	if ((em->energy() * sin(em->position().theta())) >= etmin_) isMoreEt=true;
	float e3x3 = EcalClusterTools::e3x3( *(em->seed() ), EBhit_collection, &(*topology));
	if (em->clustersSize() == 1 ) {
           const reco::BasicCluster leseed = *(em->seed());
           const std::vector< std::pair<DetId, float> > listcristaux = leseed.hitsAndFractions();
           if (nbXtals < listcristaux.size() ) {
		nbXtals = listcristaux.size();
	   }	
        }

	//if ( (e3x3/em->rawEnergy())==1 ) isGoodR9 = true;
     }
     for(reco::SuperClusterCollection::const_iterator em = superClustersEE->begin();
        em != superClustersEE->end(); ++em) {
        const std::vector< std::pair<DetId, float> > lesCristaux = em->hitsAndFractions();
        if (lesCristaux.size() == 20) isGoodSC = true;
        if ((em->energy() * sin(em->position().theta())) >= etmin_) isMoreEt=true;
        float e3x3 = EcalClusterTools::e3x3( *(em->seed() ), EEhit_collection, &(*topology))+em->preshowerEnergy();
        if ( (e3x3/em->rawEnergy())==1 ) isGoodR9 = true;
        if (em->clustersSize() == 1 ) {
           const reco::BasicCluster leseed = *(em->seed());
           const std::vector< std::pair<DetId, float> > listcristaux = leseed.hitsAndFractions();
           if (nbXtals < listcristaux.size() ) {
        //        nbXtals = listcristaux.size();
           }
        }
     }*/
     if (!(runNumber_ == iEvent.id().run())) return false; 
//   if ((keepOnly20_)&&(!(isGoodSC))) return false;
    if (!(nEvent_==-1)) {
	bool isInEvents = false;
  	int eventNumber = iEvent.id().event();
	for (int i = 0 ; i < nEvent_ ; i++){
		if (eventNumbers_[i]==eventNumber){
		 isInEvents = true;
		}
        } 
        if ( !(isInEvents) ) return false;
     }

//       if ((!(isMoreEt))||(!(isGoodR9))) return false;
//       if ((!(isMoreEt))) return false;
//       if (Xtalmin_>=nbXtals) return false;  	*/
     
    return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HughFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HughFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HughFilter);
