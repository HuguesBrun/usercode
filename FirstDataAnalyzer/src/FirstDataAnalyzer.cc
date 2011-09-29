// -*- C++ -*-
//
// Package:    FirstDataAnalyzer
// Class:      FirstDataAnalyzer
// 
/**\class FirstDataAnalyzer FirstDataAnalyzer.cc hugues/FirstDataAnalyzer/src/FirstDataAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Hugues Louis Brun
//         Created:  Mon Dec 21 16:16:59 CET 2009
// $Id: FirstDataAnalyzer.cc,v 1.2 2010/03/29 14:47:50 hbrun Exp $
//
//



//
// static data member definitions
//

//
// constructors and destructor
//


#include "../interface/FirstDataAnalyzer.h"


FirstDataAnalyzer::FirstDataAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
        isDoRecHits_ = iConfig.getParameter<bool>("readRecHits");
        isDoTTflag_  = iConfig.getParameter<bool>("readTTinfos");
	isDoHFrecHits_ = iConfig.getParameter<bool>("readHFrecHits");
	isMCTruth_ = iConfig.getParameter<bool>("readMCTruth");
	isPhotonMCTruth_ = iConfig.getParameter<bool>("readPhotonMCTruth");
	deltaRMax_ = iConfig.getParameter<double>("deltaRMax");
	triggerL1Tag_ = iConfig.getParameter<edm::InputTag>("L1triggerResults");
        HFrecHitsCollection_ = iConfig.getParameter<edm::InputTag>("HFrecHitsCollection");
	barrelEcalHits_ = iConfig.getParameter<edm::InputTag>("barrelEcalHits");
	endcapEcalHits_ = iConfig.getParameter<edm::InputTag>("endcapEcalHits");
        barrelSrpFlagsCollection_ = iConfig.getParameter<std::string>("barrelSrpFlagsCollection");
        endcapSrpFlagsCollection_ = iConfig.getParameter<std::string>("endcapSrpFlagsCollection");
        barrelClusterCollection_ = iConfig.getParameter<std::string>("barrelClusterCollection");
        barrelClusterProducer_   = iConfig.getParameter<std::string>("barrelClusterProducer");
        endcapClusterCollection_ = iConfig.getParameter<std::string>("endcapClusterCollection");
        endcapClusterProducer_   = iConfig.getParameter<std::string>("endcapClusterProducer");
	barrelCorrectedSuperClusterCollection_ = iConfig.getParameter<std::string>("barrelCorrectedSuperClusterCollection");
	barrelCorrectedSuperClusterProducer_   = iConfig.getParameter<std::string>("barrelCorrectedSuperClusterProducer");
        endcapCorrectedSuperClusterCollection_ = iConfig.getParameter<std::string>("endcapCorrectedSuperClusterCollection");
        endcapCorrectedSuperClusterProducer_   = iConfig.getParameter<std::string>("endcapCorrectedSuperClusterProducer");
	photonCollection_		       = iConfig.getParameter<std::string>("photonCollection");
	electronCollection_		       = iConfig.getParameter<std::string>("electronCollection");	
        MCParticlesCollection_                 = iConfig.getParameter<std::string>("MCParticlesCollection");       	
	Geant4SimHitsCollection_	       = iConfig.getParameter<std::string>("Geant4SimHitsCollection");
	outputFile_			       = iConfig.getParameter<std::string>("outputFile");

        rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms

	thePhotonMCTruthFinder_ = new PhotonMCTruthFinder(); // build the Nancy MC truth 

}




FirstDataAnalyzer::~FirstDataAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete rootFile_; 

}


//
// member functions 
//
void
FirstDataAnalyzer::mySuperClusterAnalyzer(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::InputTag triggerL1Tag_, edm::InputTag HFrecHitsCollection_, edm::InputTag ecalHits_, std::string srpFlagsCollection_, std::string clusterCollection_,std::string clusterProducer_,std::string correctedSuperClusterCollection_,std::string correctedSuperClusterProducer_,std::string photonCollection_,std::string electronCollection, std::string MCParticlesCollection_, bool isBarrel)
{
   using namespace edm;  // the framework classes
   using namespace std;


//getting the cms geometry and the topology of CMS and the intercalib constant
  const CaloGeometry* caloGeom;
  iSetup.get<CaloGeometryRecord>().get(theCaloGeom_);
  const CaloSubdetectorGeometry* geom=theCaloGeom_->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);//EcalBarrel = 1
  caloGeom = theCaloGeom_.product();

  iSetup.get<CaloTopologyRecord>().get(theCaloTopology_);
  const CaloTopology *topology = theCaloTopology_.product();

  edm::ESHandle<EcalIntercalibConstants> ic;
  iSetup.get<EcalIntercalibConstantsRcd>().get(ic);

  edm::ESHandle<EcalClusterLocalContCorrParameters> theParams;
  iSetup.get<EcalClusterLocalContCorrParametersRcd>().get(theParams);
  const EcalClusterLocalContCorrParameters *myParams = theParams.product();

  
  edm::ESHandle<EcalClusterCrackCorrParameters> theParamsCrack;
  iSetup.get<EcalClusterCrackCorrParametersRcd>().get(theParamsCrack);
  const EcalClusterCrackCorrParameters *myParamsCrack = theParamsCrack.product();

//  std::vector<float> myParams = theParams->params();

// Getting the RecHits collection
   edm::Handle<EcalRecHitCollection> rhcHandle_;
   iEvent.getByLabel(ecalHits_,rhcHandle_);  
   if (!(rhcHandle_.isValid())) {
     edm::LogInfo("EnergyScaleAnalyzer") << "could not get a handle on the EcalRecHitCollection!";
     return;
   }
   const EcalRecHitCollection *hit_collection = rhcHandle_.product();


// Getting the SRP flags collection
/*   edm::Handle<EBSrFlagCollection> flagHandleEB_;
   edm::Handle<EESrFlagCollection> flagHandleEE_;
   if (isBarrel) {
	iEvent.getByType(flagHandleEB_);
     //iEvent.getByLabel(srpFlagsCollection_,flagHandleEB_);
//     const EBDigiCollection *flagEB_collection = flagHandleEB_.product();
     }
     else {
     iEvent.getByLabel(srpFlagsCollection_,flagHandleEE_);
  //   const EEDigiCollection *flagEE_collection = flagHandleEE_.product();
     }
   edm::ESHandle<EcalTrigTowerConstituentsMap> hTriggerTowerMap;
   iSetup.get<IdealGeometryRecord>().get(hTriggerTowerMap);
   if (!(hTriggerTowerMap.isValid())) {
     edm::LogInfo("EnergyScaleAnalyzer") << "could not get a handle on the Trigger Tower Map collection";
   }
   const EcalTrigTowerConstituentsMap* triggerTowerMap = hTriggerTowerMap.product();*/

   

// Getting the Basic Clusters
  Handle<reco::BasicClusterCollection> pBasicClusters;
  try {
    iEvent.getByLabel(clusterProducer_, clusterCollection_, pBasicClusters);
  } catch ( cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer")
      << "L210 Error! can't get collection with label "
      << clusterCollection_.c_str() ;
  }
  //const reco::BasicClusterCollection* basicClusters = pBasicClusters.product();
// Getting the Super Clusters
  Handle<reco::SuperClusterCollection> pSuperClusters;
  try {
    iEvent.getByLabel(correctedSuperClusterProducer_,correctedSuperClusterCollection_, pSuperClusters);
  } catch ( cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer")
      << "L210 Error! can't get collection with label "
      << correctedSuperClusterCollection_.c_str() ;
  }
  const reco::SuperClusterCollection* superClusters = pSuperClusters.product();
// Getting the Photons
  Handle<reco::PhotonCollection> pPhotons;
  try {
    iEvent.getByLabel(photonCollection_, pPhotons);
  } catch (cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer")
      << "L210 Error! can't get collection with label "
      << photonCollection_.c_str() ;
  }
  const reco::PhotonCollection* photons = pPhotons.product();

// Getting the electrons

  Handle<reco::GsfElectronCollection> pElectrons;
  try {
    iEvent.getByLabel(electronCollection,pElectrons);
  } catch (cms::Exception& ex ) {
       edm::LogError("RecoPhotonEnergyScaleAnalyzer")
      << "L210 Error! can't get collection with label "
      << electronCollection_.c_str() ;
  }
  const reco::GsfElectronCollection* electrons = pElectrons.product();

// Getting the MC particles if MC datas  
  Handle<HepMCProduct> hepMC;
  if (isMCTruth_) {
     iEvent.getByLabel( MCParticlesCollection_, hepMC ) ;
     if ( !(hepMC.isValid())) {
       LogInfo("RecoPhotonEnergyScaleAnalyzer") << "Could not get MC Product!";
       return;
     }
   }
// info on the MC event
   std::vector< Handle< HepMCProduct > > evtHandles ;
     iEvent.getManyByType( evtHandles ) ;
     for ( unsigned int i=0; i<evtHandles.size(); ++i) {
     if ( evtHandles[i].isValid() ) {
        const HepMC::GenEvent* evt = evtHandles[i]->GetEvent() ;
      
      // take only 1st vertex for now - it's been tested only of PGuns...
      //
        HepMC::GenEvent::vertex_const_iterator vtx = evt->vertices_begin() ;
        if ( evtHandles[i].provenance()->moduleLabel() == MCParticlesCollection_ ) {
          //Corrdinates of Vertex w.r.o. the point (0,0,0)
          tree_.mc_vtxX = 0.1*(*vtx)->position().x();      
          tree_.mc_vtxY = 0.1*(*vtx)->position().y();
          tree_.mc_vtxZ = 0.1*(*vtx)->position().z();  
        }
      }
     }
   //std::cout << "primary vertex x = " << tree_.mc_vtxX << " y = " << tree_.mc_vtxY << " z = " << tree_.mc_vtxZ << std::endl; 
   math::XYZVector vtxPosition(tree_.mc_vtxX, tree_.mc_vtxY, tree_.mc_vtxZ);
// L1 technical trigger
  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel(triggerL1Tag_, gtReadoutRecord);
  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();


// HF collection 
   edm::Handle<HFRecHitCollection> HFrecHits;
   float eHfNeg = 0;
   float eHfPos = 0;
   float eHfNegTime = 0;
   float eHfPosTime = 0; 
   double alpha = 0;
   double alphaTime = 0; 	 
   int eHfNcounts = 0; 
   int eHfPcounts = 0;
   if (!(isBarrel)) {
       try{
          iEvent.getByLabel(HFrecHitsCollection_, HFrecHits);	
       } catch (cms::Exception& ex) {
         edm::LogWarning("RecoPhotonEnergyScaleAnalyzer")
         << "I can t read the HF rec hits collection "; 
       } 	
// calculate the Energy in the HF , the nb of count and the timming
       for(size_t ihit = 0 ; ihit < HFrecHits->size() ; ++ihit){
	  const HFRecHit h = (*HFrecHits)[ihit];
          double energy = h.energy();
	  double time  = h.time();
          const HcalDetId id(h.id());
          const CaloCellGeometry *thisHFCell = caloGeom->getGeometry(h.id());
          GlobalPoint HFposition = thisHFCell->getPosition();
	  if (isDoHFrecHits_){
	     int depthOther = 1; 
	     if (id.depth() == 1 ) depthOther = 2; //Now depthOther is the depth of the other fiber
	     HcalSubdetector subdet = id.subdet(); // recup the name of the subdetector  	
	     const HcalDetId otherId(subdet, id.ieta(), id.iphi(), depthOther); // id of the other fiber  	
	     float otherEnergy = 0;
	     float otherTime = 0;		
	     for(size_t ihit2 = 0 ; ihit2 < HFrecHits->size() ; ++ihit2){
		const HFRecHit h2 = (*HFrecHits)[ihit2];
	     	const HcalDetId id2(h2.id());	
		/*    std::cout << "coucou phi = " << otherId.iphi() << " eta = " << otherId.ieta() << " depth = " << otherId.depth() <<  std::endl;
		    std::cout << "coucou phi = " << id2.iphi() << " eta = " << id2.ieta() << " depth = " << id2.depth() <<  std::endl;
		    std::cout << " -----------------------------------------------------------------------------------------------------" << std::endl;					*/
		if (otherId == id2){ 
//		    std::cout << "coucou phi = " << id.iphi() << " eta = " << id.ieta() << " depth = " << id.depth() <<  std::endl;
//		    std::cout << "coucou phi = " << id2.iphi() << " eta = " << id2.ieta() << " depth = " << id2.depth() <<  std::endl;
		    otherEnergy = h2.energy();
		    otherTime = h2.time();
                    break; 
		}
             }
             alpha = (energy-otherEnergy)/(energy+otherEnergy);
             alphaTime = (energy*time-otherEnergy*otherTime)/(energy*time+otherEnergy*otherTime);
             if (id.depth() == 1) {alpha = -alpha; alphaTime = -alphaTime;}
//             std::cout << "alpha = " << alpha << " alphaTime = " << alphaTime << std::endl;
	     treeHF_.hf_eta = HFposition.eta();
	     treeHF_.hf_phi = HFposition.phi(); 	
	     treeHF_.hf_ieta = id.ieta();
	     treeHF_.hf_iphi = id.iphi();
	     treeHF_.hf_energy = energy;
	     treeHF_.hf_time = time;	 
	     treeHF_.hf_depth = id.depth();
             treeHF_.hf_alphaRatio = alpha;
	     treeHF_.hf_alphaRatioTimed = alphaTime;	
	     myHFTree_->Fill();	
	  }
          if (fabs(HFposition.eta())>3.5) continue;
	  if (energy < 0)  continue; 
          if (id.zside()<0) {
             eHfNeg     += energy;
             eHfNegTime += energy * time;
             if (energy>3) ++eHfNcounts;
          } else {
             eHfPos     += energy;
             eHfPosTime += energy * time;
             if (energy>3) ++eHfPcounts;
          }
       }	
  } 

// ref to the RecHits in a SC or BC
  std::vector< DetId>  recHitsInSC;
  std::vector< DetId>  recHitsInBC;
// start of the loop on the SC
  for(reco::SuperClusterCollection::const_iterator em = superClusters->begin();
        em != superClusters->end(); ++em) {
      if (isBarrel) nbSCEB_++;
	else {
	if (em->position().eta()>0) nbSCEEP_++;
	  else nbSCEEM_++;	
        }
      // Trigger L1
      if (technicalTriggerWordBeforeMask.at(0)) tree_.techTrigger0 = 1;
        else tree_.techTrigger0 = 0;
      if (technicalTriggerWordBeforeMask.at(40)) tree_.techTrigger40 = 1;
        else tree_.techTrigger40 = 0;
      if (technicalTriggerWordBeforeMask.at(41)) tree_.techTrigger41 = 1;
        else tree_.techTrigger41 = 0;
      if (technicalTriggerWordBeforeMask.at(36)) tree_.techTrigger36 = 1;
        else tree_.techTrigger36 = 0;
      if (technicalTriggerWordBeforeMask.at(37)) tree_.techTrigger37 = 1;
        else tree_.techTrigger37 = 0;
      if (technicalTriggerWordBeforeMask.at(38)) tree_.techTrigger38 = 1;
        else tree_.techTrigger38 = 0;
      if (technicalTriggerWordBeforeMask.at(39)) tree_.techTrigger39 = 1;
        else tree_.techTrigger39 = 0;

      // event caracteristics
      tree_.eventRef = iEvent.id().event();
      tree_.runNum = iEvent.id().run();
      tree_.bx = iEvent.bunchCrossing();
      tree_.orbite = iEvent.orbitNumber();
      tree_.triggerType = iEvent.experimentType();
      tree_.lumiBlock = iEvent.luminosityBlock();

     // is the SC in cracks ??
      tree_.em_isInCrack = 0;
      double emAbsEta = fabs(em->position().eta());
      // copied from RecoEgama/EgammaElectronAlgos/src/EgammaElectronClassification.cc
      if ( emAbsEta < 0.018 ||
         (emAbsEta > 0.423 && emAbsEta < 0.461) ||
         (emAbsEta > 0.770 && emAbsEta < 0.806) ||
         (emAbsEta > 1.127 && emAbsEta < 1.163) ||
         (emAbsEta > 1.460 && emAbsEta < 1.558) )
            tree_.em_isInCrack = 1;
      // part of CMS ECAL
      if (isBarrel) tree_.em_barrelOrEndcap = 1;
      else tree_.em_barrelOrEndcap = 0;	 
      // energy and position of the SC
      tree_.em_e     = em->energy();
      if (isBarrel) tree_.em_eRAW  = em->rawEnergy();
      else tree_.em_eRAW  = em->rawEnergy() + em->preshowerEnergy();
      tree_.em_et    = em->energy() * sin(em->position().theta());
      tree_.em_etRAW = tree_.em_eRAW * sin(em->position().theta());
      tree_.em_phi   = em->position().phi();
      tree_.em_eta   = em->position().eta();
      tree_.em_theta = em->position().theta();
      tree_.em_e5x5 = EcalClusterTools::e5x5( *(em->seed() ), hit_collection, &(*topology));
      tree_.em_e2x2 = EcalClusterTools::e2x2( *(em->seed() ), hit_collection, &(*topology));
      // shape of the super cluster
      tree_.em_pw1 = em->phiWidth();
      tree_.em_ew1 = em->etaWidth();
      tree_.em_sigmaetaeta = (EcalClusterTools::scLocalCovariances(*(em), hit_collection, &(*topology), 4.2)).at(0);
      tree_.em_sigmaphiphi = (EcalClusterTools::scLocalCovariances(*(em), hit_collection, &(*topology), 4.2)).at(2);
      tree_.em_BCsigmaietaieta = (EcalClusterTools::covariances(*(em->seed()), hit_collection, &(*topology), &(*caloGeom), 4.2)).at(0);
      tree_.em_BCsigmaiphiiphi = (EcalClusterTools::covariances(*(em->seed()), hit_collection, &(*topology), &(*caloGeom), 4.2)).at(2);
      tree_.em_BCsigmaetaeta = (EcalClusterTools::localCovariances(*(em->seed()), hit_collection, &(*topology), 4.2)).at(0);
      tree_.em_BCsigmaphiphi = (EcalClusterTools::localCovariances(*(em->seed()), hit_collection, &(*topology), 4.2)).at(2);
      float e3x3 = EcalClusterTools::e3x3( *(em->seed() ), hit_collection, &(*topology));
      float r9 = e3x3/em->rawEnergy();
      tree_.em_r9 = r9;
      const std::pair<DetId, float> seedRecHits = EcalClusterTools::getMaximum( (*(em->seed())).hitsAndFractions(), hit_collection);
      float eSeed = EcalClusterTools::recHitEnergy( seedRecHits.first, hit_collection);
      tree_.em_seedBCPhi = em->seed()->phi();
      tree_.em_seedBCEta = em->seed()->eta();
     // std::cout << "BC phi = " << em->seed()->phi() << " eta = " << em->seed()->eta() << std::endl;
      tree_.em_seedEnergy = eSeed;		
      float r19 = -10;
      if (!(e3x3==0)) r19=eSeed/e3x3;
      tree_.em_rookEnergy = myCaloTools::rookEnergy( hit_collection, &(*topology), seedRecHits.first);
      tree_.em_seedChi2 = myCaloTools::recupChi2(hit_collection, seedRecHits.first);
      tree_.em_seedTime = myCaloTools::recupTime(hit_collection, seedRecHits.first);
      tree_.em_seedFlag = myCaloTools::recupFlag(hit_collection, seedRecHits.first);
      tree_.em_hasBadSrpFlag = 0;
      const CaloCellGeometry *this_cell = caloGeom->getGeometry(seedRecHits.first);
      GlobalPoint position = this_cell->getPosition();
      tree_.em_seedPhi = position.phi();
      tree_.em_seedEta = position.eta();
      //std::cout << "la position du Xtal = Phi = " << position.phi() << " eta = " << position.eta() << std::endl;
      //std::cout << "la position du xtal = x  = " << position.x() << " y= " << position.y() << " z = " << position.z() << " r= " << sqrt(position.x()*position.x()+position.y()*position.y()) << std::endl;	
	math::XYZPoint coucou(1,2,3);
      if (isBarrel) {
//      const EBSrFlagCollection *flagEB_collection = flagHandleEB_.product();
      EBDetId *test = new EBDetId((seedRecHits.first));
      tree_.em_seedIphi = test->iphi();
      tree_.em_seedIeta = test->ieta();
//      std::cout << tree_.em_seedIphi << " ieta = " << tree_.em_seedIeta << std::endl;
  //    EBSrFlagCollection::const_iterator leFlag = flagEB_collection->find(triggerTowerMap->towerOf(*test));
    //  tree_.em_seedSrpFlag = leFlag->value();
 /*     for (unsigned int iteScRh = 0 ; iteScRh < em->hitsAndFractions().size() ; iteScRh++){
        EBDetId *leId = new EBDetId((em->hitsAndFractions()[iteScRh]).first);
        EBSrFlagCollection::const_iterator leFlag = flagEB_collection->find(triggerTowerMap->towerOf(*leId));
        if ((leFlag->value()==5)||(leFlag->value()==7)) {
 	   tree_.em_hasBadSrpFlag = 1;
           break;
        }
      }*/

      } 
      else {
      bool isEndcapBadSC = false;
      for (unsigned int iteScRh = 0 ; iteScRh < em->hitsAndFractions().size() ; iteScRh++){
        EEDetId *leId = new EEDetId((em->hitsAndFractions()[iteScRh]).first);
        if ((leId->ix()>40)&&(leId->iy()>40)&&(leId->ix()<61)&&(leId->iy()<61)) isEndcapBadSC = true;
      } 
      if (isEndcapBadSC) continue;
//      const EESrFlagCollection *flagEE_collection = flagHandleEE_.product();
      EEDetId *test = new EEDetId((seedRecHits.first));     
      tree_.em_seedIx = test->ix();
      tree_.em_seedIy = test->iy();
      tree_.em_seedZside = test->zside();
      const int scEdge = 5;
  /*    EcalScDetId id = EcalScDetId(((*test).ix()-1)/scEdge+1, ((*test).iy()-1)/scEdge+1, (*test).zside()); 
      EESrFlagCollection::const_iterator leFlag = flagEE_collection->find(id);
      tree_.em_seedSrpFlag = leFlag->value();*/
      }
      if (!(eSeed==0)) tree_.em_fracRook = tree_.em_rookEnergy/eSeed;
      tree_.em_r19 = r19;
      tree_.em_br1 = tree_.em_pw1/tree_.em_ew1;
      tree_.em_nBC   = em->clustersSize();
      tree_.em_nbcrystal = (em->hitsAndFractions()).size();
      const reco::BasicCluster leseed = *(em->seed());
      const std::vector< std::pair<DetId, float> > listcristaux = leseed.hitsAndFractions();
      tree_.em_nbcrystalSEED = listcristaux.size();
      if (em->clustersSize() == 1 ) {
        DetId seedXtalId = EcalClusterTools::getMaximum( leseed.hitsAndFractions(), hit_collection).first;
        tree_.em_isholeinfirst = myCaloTools::contourinf(leseed, seedXtalId, isBarrel);
        tree_.em_isholeinsecond = myCaloTools::contoursup(leseed, seedXtalId, isBarrel);
      }
      else{
        tree_.em_isholeinfirst = -10;
        tree_.em_isholeinsecond = -10;
      }
      for (unsigned int iteScRh = 0 ; iteScRh < em->hitsAndFractions().size() ; iteScRh++){
	recHitsInSC.push_back((em->hitsAndFractions()[iteScRh]).first); 
      }
     // applying the corrections
      float energyWithEtaCorrection = 0;
      if (isBarrel) energyWithEtaCorrection = tree_.em_eRAW/fcorr::f5x5((emAbsEta*(5/0.087)));
      else energyWithEtaCorrection = tree_.em_eRAW; 	
      float pTWithEtaCorrection     = energyWithEtaCorrection * sin(tree_.em_theta);
      tree_.emCorrEta_e = energyWithEtaCorrection;
      tree_.emCorrEta_et = pTWithEtaCorrection;
      if (isBarrel) tree_.em_e5x5corrEta = tree_.em_e5x5/fcorr::f5x5((TMath::Abs(em->eta()*(5/0.087))));
	else tree_.em_e5x5corrEta = tree_.em_e5x5;
      if (isBarrel) tree_.emCorrBR1_e  = fcorr::electron_br1(tree_.em_br1, tree_.emCorrEta_e);
      else tree_.emCorrBR1_e  = fcorr_ee::electron_br1(tree_.em_br1, tree_.emCorrEta_e);
      tree_.emCorrBR1_et = tree_.emCorrBR1_e * sin(em->position().theta());
      if (isBarrel) tree_.emCorrBR1Full_et = fcorr::electron_br1_complete(tree_.emCorrBR1_et, emAbsEta);
      else tree_.emCorrBR1Full_et = fcorr_ee::electron_br1_complete(tree_.emCorrBR1_et, emAbsEta);
      tree_.emCorrBR1Full_e  = tree_.emCorrBR1Full_et/sin(em->position().theta());
    // apply the Locci Local Cont corrections
      if (isBarrel) {
	 tree_.em_isEtaInner = fLocalCorr::leftOrRightEta(seedRecHits.first, hit_collection);
         tree_.em_isPhiInner = fLocalCorr::leftOrRightPhi(seedRecHits.first, hit_collection);
         tree_.em_logERatioEta = fLocalCorr::logRatioEta(seedRecHits.first, hit_collection, tree_.em_isEtaInner, &(*topology));
         tree_.em_logERatioPhi = fLocalCorr::logRatioPhi(seedRecHits.first, hit_collection, tree_.em_isPhiInner, &(*topology));
//f (tree_.em_r9>0.93) std::cout << "ieta = " << tree_.em_seedIeta << std::endl;
	 int umbParam = 84 - fabs(tree_.em_seedIeta);
	 float corrUmb = 0.992911+0.00017885*umbParam-1.04615e-06*umbParam*umbParam;
	 if (corrUmb>1) corrUmb=1;
//if (tree_.em_r9>0.93)  std::cout << "corrUmb = " << corrUmb << std::endl;
	 tree_.em_e5x5Umbrella = tree_.em_e5x5/corrUmb;
	 tree_.em_corrLocci = fLocalCorr::SLocciLocalEtaCrackCorrections(tree_.em_seedIphi, tree_.em_seedIeta, tree_.em_logERatioEta, tree_.em_isEtaInner);

	 tree_.em_corrLocalTFactorEta = fLocalCorr::theTourneurLocalCorrection( (*em), &(*geom) , myParams,1 );
	 tree_.em_corrLocalTFactorPhi = fLocalCorr::theTourneurLocalCorrection( (*em), &(*geom) , myParams,0 );
         tree_.em_corrCrackTFactorEta = fLocalCorr::theTourneurCrackCorrection( (*em), &(*geom) , myParamsCrack, 1);
         tree_.em_corrCrackTFactorPhi = fLocalCorr::theTourneurCrackCorrection( (*em), &(*geom) , myParamsCrack, 0 );

	 tree_.em_e5x5corrLocalT = tree_.em_corrLocalTFactorEta* tree_.em_corrLocalTFactorPhi * tree_.em_e5x5corrEta;
         tree_.em_e5x5corrCrackT = tree_.em_e5x5corrLocalT * tree_.em_corrCrackTFactorEta * tree_.em_corrCrackTFactorPhi;  

//	 std::cout << "log Ratio Eta = " << tree_.em_logERatioEta << " log ratio phi = " << tree_.em_logERatioPhi << std::endl;
	// test on S Touneur Corrections : 
	//EcalClusterCrackCorrection *myClass = new EcalClusterCrackCorrection(iSetup);
	//alClusterLocalContCorrection *localContCorrection = new EcalClusterLocalContCorrection(iSetup);
       //localContCorrection->EcalClusterLocalContCorrection(iSetup);        
//	std::cout << "localCorrection factor" << localContCorrection->getValue((*em), 1) << std::endl;
//        std::cout << fLocalCorr::theTourneurLocalCorrection( (*em), &(*geom) , myParams ) << std::endl;
//	std::cout << fLocalCorr::theTourneurCrackCorrection( (*em), &(*geom) , myParamsCrack ) << std::endl;
      }		 	

      // now check if the SC became a photon
      reco::PhotonCollection::const_iterator myphoton;
      int nbphoton = 0;
      for(reco::PhotonCollection::const_iterator aPho = photons->begin(); aPho != photons->end(); aPho++){
//	std::cout << "coucou" << aPho->eta() << " " << aPho->phi() << " " << aPho->energy() * sin(aPho->theta()) << std::endl;
      const reco::SuperCluster *theSC = aPho->superCluster().get();
         if ((theSC->position().eta() == em->position().eta())||(theSC->position().phi() == em->position().phi())){
              myphoton = aPho;
              nbphoton++;
         }
      }
      if (nbphoton==1){ // Yes there is a photon
         tree_.em_isPhoton=1;
         tree_.pho_e = myphoton->energy();
         tree_.pho_et = myphoton->energy() * sin(myphoton->theta());
         tree_.pho_eta = myphoton->eta();
         tree_.pho_phi = myphoton->phi();
         tree_.pho_theta = myphoton->theta();
         tree_.pho_r9 = myphoton->r9();
	 tree_.pho_isLeading = myCaloTools::isLeadingPhoton(tree_.pho_et, photons);
      	 //std::cout << tree_.pho_et << " " << tree_.pho_r9 <<" " << tree_.pho_eta << std::endl;
//²	if (isBarrel&&(tree_.pho_isLeading)&&(tree_.pho_r9>0.93)) std::cout << "coucou ya un photon !! iphi  = " << tree_.em_seedIphi << " ieta = " << tree_.em_seedIeta << "energie tranverse = " << tree_.pho_et << " " << tree_.eventRef << std::endl;

      
      }
      else if (nbphoton==0){ // no photon from that SC
         tree_.em_isPhoton=0;
      }
      else {
         cout << "strange number of photons = " << nbphoton << endl;
      }
      // is the photon converted ?

      if (nbphoton==1) {
         reco::ConversionRefVector convrefvector = myphoton->conversions();
         if (convrefvector.size() == 1){
            tree_.pho_isConverted  = 1;
            reco::ConversionRef lareference = convrefvector[0];
            const reco::Conversion *laconversion = lareference.get();
            tree_.pho_nTracks = laconversion->nTracks();
            tree_.pho_EoverP = laconversion->EoverP();
            tree_.pho_Rconv = sqrt(laconversion->conversionVertex().x()*laconversion->conversionVertex().x()+laconversion->conversionVertex().y()*laconversion->conversionVertex().y());
            tree_.pho_Zconv = laconversion->conversionVertex().z();
            tree_.pho_Xconv = laconversion->conversionVertex().x();
            tree_.pho_Yconv = laconversion->conversionVertex().y();
         }
         else {
            tree_.pho_isConverted  = 0;
         }
      }
      // now check is it is an electron 
      reco::GsfElectronCollection::const_iterator myElectron;
      int nbelectron = 0;
      for(reco::GsfElectronCollection::const_iterator aElec = electrons->begin();  aElec!=electrons->end(); aElec++){
	const reco::SuperCluster *theElecSC = aElec->superCluster().get();
         if ((theElecSC->position().eta() == em->position().eta())||(theElecSC->position().phi() == em->position().phi())){
              myElectron = aElec;
              nbelectron++;
         }
      }
      if (nbelectron == 1) { //Yes, there is an electron
	tree_.em_isElectron = 1;
	tree_.ele_e = myElectron->energy();
	tree_.ele_et = myElectron->energy() * sin(myElectron->theta());
	tree_.ele_phi = myElectron->phi();
	tree_.ele_eta = myElectron->eta();
	tree_.ele_theta = myElectron->theta();
	tree_.ele_charge = myElectron->charge();	
	tree_.ele_mass = myElectron->mass();
	tree_.ele_EoverP = myElectron->eSuperClusterOverP();
      }
      else if (nbelectron == 0 ) {
	tree_.em_isElectron = 0;
      }		
      else {
	tree_.em_isElectron = 2;
	std::cout << "strange number of electron " << nbelectron << std::endl;
      }	
      // now check the MC truth !!!!!
      if (isMCTruth_) {
	 const HepMC::GenEvent* genEvent = hepMC->GetEvent();
         HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
         std::vector<HepMC::GenParticle* > mcParticles;
	 for ( HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
         p != genEvent->particles_end();
         ++p ) {
               	if ( (*p)->pdg_id() != 22 && abs((*p)->pdg_id()) != 11 ) continue;
		if ( (*p)->momentum().perp() < 0.8 ) continue;
		float dR = kinem::delta_R((*p)->momentum().eta(), (*p)->momentum().phi(), em->position().eta(), em->position().phi());
		if (dR < deltaRMax_) mcParticles.push_back(*p);
	 }
	 HepMC::GenParticle* mc1 = new HepMC::GenParticle(); // on cree de tout facon MC1
	 if (mcParticles.size()==0) tree_.em_isMatchWithMC = 0;
	 else if (mcParticles.size()==1) {
	     tree_.em_isMatchWithMC = 1;
	     tree_.mc_nbMatch = 1;
	     tree_.mc_PDGType = mcParticles[0]->pdg_id();
	     tree_.mc_e = mcParticles[0]->momentum().e();
	     tree_.mc_et = mcParticles[0]->momentum().perp();
	     tree_.mc_eta = mcParticles[0]->momentum().eta();
	//std::cout << "the eta of the MC = " << tree_.mc_eta << " le phi = " << tree_.mc_phi << std::endl;
	//std::cout << "x = " << mcParticles[0]->momentum().x() << " y = " << mcParticles[0]->momentum().y() << " z= " << mcParticles[0]->momentum().z()<< std::endl;
//     float dist = 129.0 - sqrt(vtxPosition.x()*vtxPosition.x()+vtxPosition.y()*vtxPosition.y());
//		std::cout << "la distance = " << dist << std::endl;
//		if ((tree_.em_r9>0.93)&&(tree_.em_isPhoton==1)&&(fabs(tree_.em_eta)<1.45)) std::cout << "local = " <<  tree_.em_e5x5corrEta << " " << tree_.em_corrLocalTFactor << " crack = " << tree_.em_corrCrackTFactor << "e mc = " << tree_.mc_e << "ratio = " << tree_.em_e5x5corrCrackT/tree_.mc_e << std::endl;
	     float normTransverse = sqrt(mcParticles[0]->momentum().x()*mcParticles[0]->momentum().x()+mcParticles[0]->momentum().y()*mcParticles[0]->momentum().y());
             float coeffNorm = 129.0/normTransverse;
	     math::XYZVector MCdirection(coeffNorm*mcParticles[0]->momentum().x(),coeffNorm*mcParticles[0]->momentum().y(),coeffNorm*mcParticles[0]->momentum().z());
//		std::cout << "verif normalisation = " << sqrt(MCdirection.x()*MCdirection.x()+MCdirection.y()*MCdirection.y()) << std::endl;
	     math::XYZVector dirPhotonNominal = vtxPosition + MCdirection;
//	     std::cout << "long transverse = " << sqrt(dirPhotonNominal.x()*dirPhotonNominal.x()+dirPhotonNominal.y()*dirPhotonNominal.y()) << std::endl;
  // 	     std::cout << "le nouvel eta " << dirPhotonNominal.eta() << std::endl;
	     tree_.mc_etaCorrVtx = dirPhotonNominal.eta();
	     tree_.mc_phi = mcParticles[0]->momentum().phi();
	     tree_.mc_theta = mcParticles[0]->momentum().theta();
             mc1 = mcParticles[0]; 
	 }
	 else {
		tree_.mc_nbMatch = mcParticles.size();
		tree_.em_isMatchWithMC = 1;
//		HepMC::GenParticle* mc1 = mcParticles[0];
		float theDeltaRmin  = 1;
		for (unsigned int j = 1 ; j < mcParticles.size() ; j++){
			float theDeltaR = kinem::delta_R(mcParticles[j]->momentum().eta(), mcParticles[j]->momentum().phi(), em->position().eta(), em->position().phi());
			if (theDeltaR < theDeltaRmin) {
				theDeltaRmin = theDeltaR;
				mc1 = mcParticles[j];	
			}
		}
		tree_.mc_PDGType = mcParticles[0]->pdg_id();
		tree_.mc_e = mc1->momentum().e();
		tree_.mc_et = mc1->momentum().perp();
                tree_.mc_eta = mc1->momentum().eta();
                tree_.mc_phi = mc1->momentum().phi();
                tree_.mc_theta = mc1->momentum().theta();
	      }
	if ((tree_.mc_PDGType==22)&&(mcParticles.size()>0)) {
		tree_.mc_isPhoton = 1;
		if (isPhotonMCTruth_){ // maybe it's time to get the conv photons MC truth 
			std::vector<SimTrack> theSimTracks;
			std::vector<SimVertex> theSimVertices;
  
			edm::Handle<SimTrackContainer> SimTk;
			edm::Handle<SimVertexContainer> SimVtx;
			iEvent.getByLabel(Geant4SimHitsCollection_, SimTk);
			iEvent.getByLabel(Geant4SimHitsCollection_, SimVtx);

        	        theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
			theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
			std::vector<PhotonMCTruth> mcPhotons = thePhotonMCTruthFinder_->find (theSimTracks,  theSimVertices); 
			int n_conv = 0;
			tree_.mc_isConverted = 0; 
			for ( std::vector<PhotonMCTruth>::const_iterator mcPho = mcPhotons.begin(); mcPho != mcPhotons.end(); mcPho++) {
			        // is it matched?
				double eta = (*mcPho).fourMomentum().eta();
				double phi = (*mcPho).fourMomentum().phi();
  				double dR = kinem::delta_R(eta, phi, mc1->momentum().eta(), mc1->momentum().phi());
		        	//std::cout << "dR = " << dR << std::endl;
				if ( dR < 1e-6 ) {
					tree_.mc_isConverted = (*mcPho).isAConversion();
					tree_.mc_convEt = (*mcPho).fourMomentum().et();
					tree_.mc_convR = (*mcPho).vertex().perp();
					tree_.mc_convZ = (*mcPho).vertex().mag() * (*mcPho).vertex().cosTheta();
					tree_.mc_convX = (*mcPho).vertex().x();
					tree_.mc_convY = (*mcPho).vertex().y();
		/*			std::cout << "primary vertex  x = " << (*mcPho).primaryVertex().x() << " y = " << (*mcPho).primaryVertex().y() << " z = " << (*mcPho).primaryVertex().mag()*(*mcPho).primaryVertex().cosTheta() << std::endl;
					std::cout << "mother vertex  x = " << (*mcPho).motherVtx().x() << " y = " << (*mcPho).motherVtx().y() << " z = " << (*mcPho).motherVtx().mag()*(*mcPho).motherVtx().cosTheta() << std::endl;*/
				}
			        if ( dR < 0.1 ) n_conv += 1;

			}	
			tree_.mc_Nconv = n_conv;
		}
	}
      }	
      myTree_->Fill();
  }  
  const reco::BasicClusterCollection* basicClusters = pBasicClusters.product();
  for (reco::BasicClusterCollection::const_iterator bc = basicClusters->begin();
        bc != basicClusters->end(); ++bc) {
     if (isBarrel) nbBCEB_++;
        else {
		if (bc->position().eta()>0) nbBCEEP_++;
		else nbBCEEM_++;
        }
            // Trigger L1
      if (technicalTriggerWordBeforeMask.at(0)) treeBC_.techTrigger0 = 1;
        else treeBC_.techTrigger0 = 0;
      if (technicalTriggerWordBeforeMask.at(40)) treeBC_.techTrigger40 = 1;
        else treeBC_.techTrigger40 = 0;
      if (technicalTriggerWordBeforeMask.at(41)) treeBC_.techTrigger41 = 1;
        else treeBC_.techTrigger41 = 0;
      if (technicalTriggerWordBeforeMask.at(36)) treeBC_.techTrigger36 = 1;
        else treeBC_.techTrigger36 = 0;
      if (technicalTriggerWordBeforeMask.at(37)) treeBC_.techTrigger37 = 1;
        else treeBC_.techTrigger37 = 0;
      if (technicalTriggerWordBeforeMask.at(38)) treeBC_.techTrigger38 = 1;
        else treeBC_.techTrigger38 = 0;
      if (technicalTriggerWordBeforeMask.at(39)) treeBC_.techTrigger39 = 1;
        else treeBC_.techTrigger39 = 0;

      // event caracteristics
      treeBC_.eventRef = iEvent.id().event();
      treeBC_.runNum = iEvent.id().run();
      treeBC_.bx = iEvent.bunchCrossing();
      treeBC_.orbite = iEvent.orbitNumber();
      treeBC_.triggerType = iEvent.experimentType();
      treeBC_.lumiBlock = iEvent.luminosityBlock();


     // is the SC in cracks ??
      treeBC_.bc_isInCrack = 0;
      double emAbsEta = fabs(bc->position().eta());
      // copied from RecoEgama/EgammaElectronAlgos/src/EgammaElectronClassification.cc
      if ( emAbsEta < 0.018 ||
         (emAbsEta > 0.423 && emAbsEta < 0.461) ||
         (emAbsEta > 0.770 && emAbsEta < 0.806) ||
         (emAbsEta > 1.127 && emAbsEta < 1.163) ||
         (emAbsEta > 1.460 && emAbsEta < 1.558) )
            treeBC_.bc_isInCrack = 1;
     // part of CMS ECAL
      if (isBarrel) treeBC_.bc_barrelOrEndcap = 1;
      else treeBC_.bc_barrelOrEndcap = 0;
     //infos on the BC
      treeBC_.bc_e     = bc->energy();	
      treeBC_.bc_et    = bc->energy() * sin(bc->position().theta());
      treeBC_.bc_phi   = bc->position().phi();
      treeBC_.bc_eta   = bc->position().eta();
      treeBC_.bc_r9    = EcalClusterTools::e5x5( *(bc), hit_collection, &(*topology))/bc->energy();  
      const std::vector< std::pair<DetId, float> > listcristaux = bc->hitsAndFractions();
      treeBC_.bc_nbcrystal = listcristaux.size();
      for (unsigned int iteBcRh = 0 ; iteBcRh < listcristaux.size() ; iteBcRh++){
        recHitsInBC.push_back((bc->hitsAndFractions()[iteBcRh]).first);
      }
      // check if the BC is in a SC	
      int refSC = 0; 
      for(reco::SuperClusterCollection::const_iterator em = superClusters->begin();
        em != superClusters->end(); ++em) {
	refSC++;
	treeBC_.bc_isInSC = 0;
	treeBC_.bc_refSC = 0;
	for (reco::CaloCluster_iterator embc = em->clustersBegin() ; embc != em->clustersEnd() ; embc++){
             if (embc->get()->position().phi() == bc->position().phi() ) {
		treeBC_.bc_isInSC = 1;
                treeBC_.bc_refSC = refSC;
	     }
        }	
      }
      myBCTree_->Fill(); 
  }
    if (isDoRecHits_){
      EcalRecHitCollection::const_iterator it;
      EcalIntercalibConstantMap::const_iterator itIC;
      for (it = hit_collection->begin(); it != hit_collection->end(); it++){
           // Trigger L1
          if (technicalTriggerWordBeforeMask.at(0)) treeRecHits_.techTrigger0 = 1;
             else treeRecHits_.techTrigger0 = 0;
          if (technicalTriggerWordBeforeMask.at(40)) treeRecHits_.techTrigger40 = 1;
             else treeRecHits_.techTrigger40 = 0;
          if (technicalTriggerWordBeforeMask.at(41)) treeRecHits_.techTrigger41 = 1;
             else treeRecHits_.techTrigger41 = 0;
          if (technicalTriggerWordBeforeMask.at(36)) treeRecHits_.techTrigger36 = 1;
             else treeRecHits_.techTrigger36 = 0;
          if (technicalTriggerWordBeforeMask.at(37)) treeRecHits_.techTrigger37 = 1;
             else treeRecHits_.techTrigger37 = 0;
          if (technicalTriggerWordBeforeMask.at(38)) treeRecHits_.techTrigger38 = 1;
             else treeRecHits_.techTrigger38 = 0;
          if (technicalTriggerWordBeforeMask.at(39)) treeRecHits_.techTrigger39 = 1;
             else treeRecHits_.techTrigger39 = 0;
          if (technicalTriggerWordBeforeMask.at(42)) treeRecHits_.techTrigger42 = 1;
             else treeRecHits_.techTrigger42 = 0;
          if (technicalTriggerWordBeforeMask.at(43)) treeRecHits_.techTrigger43 = 1;
             else treeRecHits_.techTrigger43 = 0;

      // event caracteristics
         treeRecHits_.eventRef = iEvent.id().event();
         treeRecHits_.runNum = iEvent.id().run();
         treeRecHits_.bx = iEvent.bunchCrossing();
         treeRecHits_.orbite = iEvent.orbitNumber();
         treeRecHits_.triggerType = iEvent.experimentType();
         treeRecHits_.lumiBlock = iEvent.luminosityBlock();
         // calc the phi and the eta of the crystal
         const CaloCellGeometry *thisCell = caloGeom->getGeometry(it->id());
         GlobalPoint position = thisCell->getPosition();
         treeRecHits_.rh_phi = position.phi();
         treeRecHits_.rh_eta = position.eta();
         if (isBarrel) {
         treeRecHits_.rh_barrelOrEndcap = 1; nbRH_EB_++;
         const EBDetId *theID = new EBDetId(it->id());
         treeRecHits_.rh_iPhi = theID->iphi();
         treeRecHits_.rh_iEta = theID->ieta();
         treeRecHits_.rh_iSm  = theID->ism();
         treeRecHits_.rh_iSmPhi = theID->iphiSM();
         treeRecHits_.rh_iSmEta = theID->ietaSM();
         treeRecHits_.rh_iX = 0;
         treeRecHits_.rh_iY = 0;
         itIC = ic->find((*theID));
         treeRecHits_.rh_uncalibEnergy = it->energy()/(*itIC);
//         const EBSrFlagCollection *flagEB_collection = flagHandleEB_.product();
 //        EBSrFlagCollection::const_iterator leFlag = flagEB_collection->find(triggerTowerMap->towerOf(*theID));
 //        treeRecHits_.rh_srpFlag = leFlag->value();
         }
         else {
         treeRecHits_.rh_barrelOrEndcap = 0; nbRH_EE_++;
         const EEDetId *theID = new EEDetId(it->id());
         treeRecHits_.rh_iPhi = 0;
         treeRecHits_.rh_iEta = 0;
         treeRecHits_.rh_iX = theID->ix();
         treeRecHits_.rh_iY = theID->iy();
         treeRecHits_.rh_zSide = theID->zside();
         itIC = ic->find((*theID));
         treeRecHits_.rh_uncalibEnergy = it->energy()/(*itIC);
//         const EESrFlagCollection *flagEE_collection = flagHandleEE_.product();
  //       const int scEdge = 5;
//         EcalScDetId id = EcalScDetId(((*theID).ix()-1)/scEdge+1, ((*theID).iy()-1)/scEdge+1, (*theID).zside());
//         EESrFlagCollection::const_iterator leFlag = flagEE_collection->find(id);
//         treeRecHits_.rh_srpFlag = leFlag->value();
         }
         treeRecHits_.rh_energy = it->energy();
         treeRecHits_.rh_chi2 = it->chi2();
         treeRecHits_.rh_time = it->time();
         treeRecHits_.rh_flag = it->flags();
         treeRecHits_.rh_isCluster = 0;
         for (unsigned int it_bc = 0 ; it_bc < recHitsInBC.size() ; it_bc++){
            if (it->id() == recHitsInBC[it_bc] ) {
               treeRecHits_.rh_isCluster = 1;
               break;
            }
         }
         treeRecHits_.rh_isSuperCluster = 0;
         for (unsigned int it_sc = 0 ; it_sc < recHitsInSC.size() ; it_sc++){
            if (it->id() == recHitsInSC[it_sc] ) {
               treeRecHits_.rh_isSuperCluster = 1;
               break;
            }
         }
         treeRecHits_.rh_eHfNeg = eHfNeg;
	 treeRecHits_.rh_eHfPos = eHfPos;
	 treeRecHits_.rh_eHfNegTime = eHfNegTime;
	 treeRecHits_.rh_eHfPosTime = eHfPosTime;
	 treeRecHits_.rh_eHfNcounts = eHfNcounts;
	 treeRecHits_.rh_eHfPcounts = eHfPcounts;
	 myRecHitsTree_->Fill();
      } 
    }
/*   if (isBarrel){
     const EBSrFlagCollection *flagEB_collection = flagHandleEB_.product();
    EBSrFlagCollection::const_iterator leFlag;
      for ( leFlag = flagEB_collection->begin() ; leFlag != flagEB_collection->end() ; leFlag++){
	if ( fabs(leFlag->id().ieta()) < 18 ){
             if (isDoTTflag_){
                   // Trigger L1
                   if (technicalTriggerWordBeforeMask.at(0)) treeTrigger_.techTrigger0 = 1;
                       else treeTrigger_.techTrigger0 = 0;
                   if (technicalTriggerWordBeforeMask.at(40)) treeTrigger_.techTrigger40 = 1;
                       else treeTrigger_.techTrigger40 = 0;
                   if (technicalTriggerWordBeforeMask.at(41)) treeTrigger_.techTrigger41 = 1;
                       else treeTrigger_.techTrigger41 = 0;
                   if (technicalTriggerWordBeforeMask.at(36)) treeTrigger_.techTrigger36 = 1;
                       else treeTrigger_.techTrigger36 = 0;
                   if (technicalTriggerWordBeforeMask.at(37)) treeTrigger_.techTrigger37 = 1;
                       else treeTrigger_.techTrigger37 = 0;
                   if (technicalTriggerWordBeforeMask.at(38)) treeTrigger_.techTrigger38 = 1;
                       else treeTrigger_.techTrigger38 = 0;
                   if (technicalTriggerWordBeforeMask.at(39)) treeTrigger_.techTrigger39 = 1;
                       else treeTrigger_.techTrigger39 = 0;

                      // event caracteristics
                   treeTrigger_.eventRef = iEvent.id().event();
                   treeTrigger_.runNum = iEvent.id().run();
                   treeTrigger_.bx = iEvent.bunchCrossing();
                   treeTrigger_.orbite = iEvent.orbitNumber();
                   treeTrigger_.triggerType = iEvent.experimentType();
                   treeTrigger_.lumiBlock = iEvent.luminosityBlock();

		   treeTrigger_.TT_barrelOrEndcap = 1;
                   treeTrigger_.TT_iphi = leFlag->id().iphi();
                   treeTrigger_.TT_ieta = leFlag->id().ieta();
                   treeTrigger_.TT_zside = leFlag->id().zside();
                   treeTrigger_.TT_flag  = leFlag->value();
                   myTrigTree_->Fill();
                 }
		switch (leFlag->value()){
		    case 1 : nbZsTTEB_++; break;
		    case 3 : nbSrpTTEB_++; break;
                    case 0 : nbFlag0TTEB_++; break;
                    case 2 : nbFlag2TTEB_++; break;
                    case 4 : nbFlag4TTEB_++; break;
                    case 5 : nbFlag5TTEB_++; break;
                    case 6 : nbFlag6TTEB_++; break;
                    case 7 : nbFlag7TTEB_++; break;	
                }
        }
     }
   }
   else {
     const EESrFlagCollection *flagEE_collection = flagHandleEE_.product();
     EESrFlagCollection::const_iterator leFlag;
      for ( leFlag = flagEE_collection->begin() ; leFlag != flagEE_collection->end() ; leFlag++){
         if (isDoTTflag_){
                   // Trigger L1
               if (technicalTriggerWordBeforeMask.at(0)) treeTrigger_.techTrigger0 = 1;
                   else treeTrigger_.techTrigger0 = 0;
               if (technicalTriggerWordBeforeMask.at(40)) treeTrigger_.techTrigger40 = 1;
                   else treeTrigger_.techTrigger40 = 0;
               if (technicalTriggerWordBeforeMask.at(41)) treeTrigger_.techTrigger41 = 1;
                   else treeTrigger_.techTrigger41 = 0;
               if (technicalTriggerWordBeforeMask.at(36)) treeTrigger_.techTrigger36 = 1;
                   else treeTrigger_.techTrigger36 = 0;
               if (technicalTriggerWordBeforeMask.at(37)) treeTrigger_.techTrigger37 = 1;
                   else treeTrigger_.techTrigger37 = 0;
               if (technicalTriggerWordBeforeMask.at(38)) treeTrigger_.techTrigger38 = 1;
                   else treeTrigger_.techTrigger38 = 0;
               if (technicalTriggerWordBeforeMask.at(39)) treeTrigger_.techTrigger39 = 1;
                   else treeTrigger_.techTrigger39 = 0;

                      // event caracteristics
               treeTrigger_.eventRef = iEvent.id().event();
               treeTrigger_.runNum = iEvent.id().run();
               treeTrigger_.bx = iEvent.bunchCrossing();
               treeTrigger_.orbite = iEvent.orbitNumber();
               treeTrigger_.triggerType = iEvent.experimentType();
               treeTrigger_.lumiBlock = iEvent.luminosityBlock();

	       treeTrigger_.TT_barrelOrEndcap = 0;	
               treeTrigger_.TT_ix = leFlag->id().ix();
               treeTrigger_.TT_iy = leFlag->id().iy();
               treeTrigger_.TT_zside = leFlag->id().zside();
               treeTrigger_.TT_flag  = leFlag->value();
               myTrigTree_->Fill();
               }

        if ((leFlag->id().zside()==-1)){
                switch (leFlag->value()){
                    case 1 : nbZsTTEEM_++; break;
                    case 3 : nbSrpTTEEM_++; break;
                    case 0 : nbFlag0TTEEM_++; break;
                    case 2 : nbFlag2TTEEM_++; break;
                    case 4 : nbFlag4TTEEM_++; break;
                    case 5 : nbFlag5TTEEM_++; break;
                    case 6 : nbFlag6TTEEM_++; break;
                    case 7 : nbFlag7TTEEM_++; break;
                }
        }
        if ((leFlag->id().zside()==1)){
                switch (leFlag->value()){
                    case 1 : nbZsTTEEP_++; break;
                    case 3 : nbSrpTTEEP_++; break;
                    case 0 : nbFlag0TTEEP_++; break;
                    case 2 : nbFlag2TTEEP_++; break;
                    case 4 : nbFlag4TTEEP_++; break;
                    case 5 : nbFlag5TTEEP_++; break;
                    case 6 : nbFlag6TTEEP_++; break;
                    case 7 : nbFlag7TTEEP_++; break;
                }
        }


     }

   }*/
}

// ------------ method called to for each event  ------------
void
FirstDataAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;  // the framework classes
   using namespace std;
   
   nbSCEB_=0;
   nbSCEEP_=0;
   nbSCEEM_=0;
   nbBCEB_=0;
   nbBCEEP_=0;
   nbBCEEM_=0;

   nbSrpTTEB_ = 0;
   nbZsTTEB_ = 0;
   nbFlag0TTEB_ = 0;
   nbFlag2TTEB_ = 0;
   nbFlag4TTEB_ = 0;
   nbFlag5TTEB_ = 0;
   nbFlag6TTEB_ = 0;
   nbFlag7TTEB_ = 0;
   nbSrpTTEEM_ = 0;
   nbZsTTEEM_ = 0;
   nbFlag0TTEEM_ = 0;
   nbFlag2TTEEM_ = 0;
   nbFlag4TTEEM_ = 0;
   nbFlag5TTEEM_ = 0;
   nbFlag6TTEEM_ = 0;
   nbFlag7TTEEM_ = 0;
   nbSrpTTEEP_ = 0;
   nbZsTTEEP_ = 0;
   nbFlag0TTEEP_ = 0;
   nbFlag2TTEEP_ = 0;
   nbFlag4TTEEP_ = 0;
   nbFlag5TTEEP_ = 0;
   nbFlag6TTEEP_ = 0;
   nbFlag7TTEEP_ = 0;
   nbRH_EB_ = 0;
   nbRH_EE_ = 0;
   
   mySuperClusterAnalyzer(iEvent, iSetup, triggerL1Tag_, HFrecHitsCollection_, barrelEcalHits_, barrelSrpFlagsCollection_, barrelClusterCollection_, barrelClusterProducer_, barrelCorrectedSuperClusterCollection_, barrelCorrectedSuperClusterProducer_, photonCollection_, electronCollection_, MCParticlesCollection_, true);

   mySuperClusterAnalyzer(iEvent, iSetup, triggerL1Tag_, HFrecHitsCollection_, endcapEcalHits_, endcapSrpFlagsCollection_, endcapClusterCollection_, endcapClusterProducer_, endcapCorrectedSuperClusterCollection_, endcapCorrectedSuperClusterProducer_, photonCollection_, electronCollection_, MCParticlesCollection_, false );
	
// L1 technical trigger
  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel(triggerL1Tag_, gtReadoutRecord);
  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();

     treeEvent_.eventRef = iEvent.id().event();
     treeEvent_.runNum = iEvent.id().run();
     treeEvent_.bx = iEvent.bunchCrossing();
     treeEvent_.orbite = iEvent.orbitNumber();
     treeEvent_.triggerType = iEvent.experimentType();
     treeEvent_.lumiBlock = iEvent.luminosityBlock();


// start of the loop on the SC
      // Trigger L1
      if (technicalTriggerWordBeforeMask.at(0)) treeEvent_.techTrigger0 = 1;
        else treeEvent_.techTrigger0 = 0;
      if (technicalTriggerWordBeforeMask.at(40)) treeEvent_.techTrigger40 = 1;
        else treeEvent_.techTrigger40 = 0;
      if (technicalTriggerWordBeforeMask.at(41)) treeEvent_.techTrigger41 = 1;
        else treeEvent_.techTrigger41 = 0;
      if (technicalTriggerWordBeforeMask.at(36)) treeEvent_.techTrigger36 = 1;
        else treeEvent_.techTrigger36 = 0;
      if (technicalTriggerWordBeforeMask.at(37)) treeEvent_.techTrigger37 = 1;
        else treeEvent_.techTrigger37 = 0;
      if (technicalTriggerWordBeforeMask.at(38)) treeEvent_.techTrigger38 = 1;
        else treeEvent_.techTrigger38 = 0;
      if (technicalTriggerWordBeforeMask.at(39)) treeEvent_.techTrigger39 = 1;
        else treeEvent_.techTrigger39 = 0;
      treeEvent_.nbSuperClusterBarrel = nbSCEB_;
      treeEvent_.nbSuperClusterEndcapM = nbSCEEM_;
      treeEvent_.nbSuperClusterEndcapP = nbSCEEP_;
      treeEvent_.nbBasicClusterBarrel = nbBCEB_;
      treeEvent_.nbBasicClusterEndcapM = nbBCEEM_;		
      treeEvent_.nbBasicClusterEndcapP = nbBCEEP_;

      treeEvent_.nbSrpTTEB = nbSrpTTEB_;
      treeEvent_.nbZsTTEB = nbZsTTEB_;
      treeEvent_.nbFlag0TTEB = nbFlag0TTEB_;
      treeEvent_.nbFlag2TTEB = nbFlag2TTEB_;
      treeEvent_.nbFlag4TTEB = nbFlag4TTEB_;
      treeEvent_.nbFlag5TTEB = nbFlag5TTEB_;
      treeEvent_.nbFlag6TTEB = nbFlag6TTEB_;
      treeEvent_.nbFlag7TTEB = nbFlag7TTEB_;
      treeEvent_.nbSrpTTEEM = nbSrpTTEEM_;
      treeEvent_.nbZsTTEEM = nbZsTTEEM_;
      treeEvent_.nbFlag0TTEEM = nbFlag0TTEEM_;
      treeEvent_.nbFlag2TTEEM = nbFlag2TTEEM_;
      treeEvent_.nbFlag4TTEEM = nbFlag4TTEEM_;
      treeEvent_.nbFlag5TTEEM = nbFlag5TTEEM_;
      treeEvent_.nbFlag6TTEEM = nbFlag6TTEEM_;
      treeEvent_.nbFlag7TTEEM = nbFlag7TTEEM_;
      treeEvent_.nbSrpTTEEP = nbSrpTTEEP_;
      treeEvent_.nbZsTTEEP = nbZsTTEEP_;
      treeEvent_.nbFlag0TTEEP = nbFlag0TTEEP_;
      treeEvent_.nbFlag2TTEEP = nbFlag2TTEEP_;
      treeEvent_.nbFlag4TTEEP = nbFlag4TTEEP_;
      treeEvent_.nbFlag5TTEEP = nbFlag5TTEEP_;
      treeEvent_.nbFlag6TTEEP = nbFlag6TTEEP_;
      treeEvent_.nbFlag7TTEEP = nbFlag7TTEEP_;
      treeEvent_.nbRH_EB = nbRH_EB_;
      treeEvent_.nbRH_EE = nbRH_EE_;

     myEventTree_->Fill();    
}


// ------------ method called once each job just before starting event loop  ------------
void 
FirstDataAnalyzer::beginJob()
{
  myTree_ = new TTree("energyScale","");
  TString treeVariables = "eventRef/I:runNum/I:bx/I:orbite/I:triggerType/I:lumiBlock/I:"; //event references
  treeVariables += "techTrigger0/I:techTrigger40/I:techTrigger41/I:techTrigger36/I:techTrigger37/I:techTrigger38/I:techTrigger39/I:"; // technical trigger 
  treeVariables += "em_isInCrack/I:em_barrelOrEndcap/I:em_e/F:em_eRAW/F:em_et/F:em_etRAW/F:em_phi/F:em_eta/F:em_theta/F:em_e5x5/F:em_e2x2/F:"; // SC infos
  treeVariables += "em_pw1/F:em_ew1/F:em_sigmaetaeta/F:em_sigmaietaieta/F:em_sigmaphiphi/F:em_sigmaiphiiphi/F:em_BCsigmaetaeta/F:em_BCsigmaietaieta/F:em_BCsigmaphiphi/F:em_BCsigmaiphiiphi/F:em_r9/F:em_r19/F:em_br1/F:em_nBC/I:em_nbcrystal/I:em_nbcrystalSEED/I:em_isholeinfirst/I:em_isholeinsecond/I:em_rookEnergy/F:em_fracRook/F:em_hasBadSrpFlag/I:em_seedBCPhi/F:em_seedBCEta/F:em_seedEnergy/F:em_seedChi2/F:em_seedTime/F:em_seedFlag/I:em_seedSrpFlag/I:em_seedPhi/F:em_seedEta/F:em_seedIphi/I:em_seedIeta/I:em_seedIx/I:em_seedIy/I:em_seedZside/I:"; // SC shape
  treeVariables += "emCorrEta_e/F:emCorrEta_et/F:em_e5x5corrEta/F:emCorrBR1_e/F:emCorrBR1_et/F:emCorrBR1Full_e/F:emCorrBR1Full_et/F:"; // SC corrections
  treeVariables += "em_isEtaInner/I:em_isPhiInner/I:em_logERatioEta/F:em_logERatioPhi/F:em_e5x5Umbrella/F:em_corrLocci/F:"; //E locci corrections
  treeVariables += "em_corrLocalTFactorEta/F:em_corrLocalTFactorPhi/F:em_corrCrackTFactorEta/F:em_corrCrackTFactorPhi/F:em_e5x5corrLocalT/F:em_e5x5corrCrackT/F:"; //S tourneur corrections	
  treeVariables += "em_isPhoton/I:pho_e/F:pho_et/F:pho_phi/F:pho_eta/F:pho_theta/F:pho_r9/F:pho_isLeading/I:pho_isConverted/I:"; // photon information
  treeVariables += "pho_nTracks/I:pho_EoverP/F:pho_Rconv/F:pho_Zconv/F:pho_Xconv/F:pho_Yconv/F:"; // converted photon
  treeVariables += "em_isElectron/I:ele_e/F:ele_et/F:ele_phi/F:ele_eta/F:ele_theta/F:ele_charge/F:ele_mass/F:ele_EoverP/F:"; //electron information
  treeVariables += "mc_vtxX/F:mc_vtxY/F:mc_vtxZ/F:em_isMatchWithMC/I:mc_PDGType/I:mc_e/F:mc_et/F:mc_eta/F:mc_phi/F:mc_theta/F:mc_nbMatch/I:mc_etaCorrVtx/F:";//general MC infos
  treeVariables += "mc_isPhoton/I:mc_isConverted/I:mc_convEt/F:mc_convR/F:mc_convR/F:mc_convX/F:mc_convY/F:mc_convZ/F:mc_Nconv/I";// MC infos for MC photons
  myTree_->Branch("energyScale",&(tree_.eventRef),treeVariables);
  

  myBCTree_ = new TTree("basicClusterTree","");
  TString treeVariables2 = "eventRef/I:runNum/I:bx/I:orbite/I:triggerType/I:lumiBlock/I:"; //event references
  treeVariables2 += "techTrigger0/I:techTrigger40/I:techTrigger41/I:techTrigger36/I:techTrigger37/I:techTrigger38/I:techTrigger39/I:"; // technical trigger
  treeVariables2 += "bc_isInCrack/I:bc_barrelOrEndcap/I:bc_e/F:bc_et/F:bc_phi/F:bc_eta/F:"; // BC infos
  treeVariables2 += "bc_isInSC/I:bc_refSC/I:bc_fraction/F:";  // ref to the SC
  treeVariables2 += "bc_r9/F:bc_nbcrystal/I"; 
  myBCTree_->Branch("basicClusterTree",&(treeBC_.eventRef),treeVariables2);

  myEventTree_ = new TTree("eventTree","");
 TString treeVariables3 = "eventRef/I:runNum/I:bx/I:orbite/I:triggerType/I:lumiBlock/I:"; //event references
  treeVariables3 += "techTrigger0/I:techTrigger40/I:techTrigger41/I:techTrigger36/I:techTrigger37/I:techTrigger38/I:techTrigger39/I:"; //trigger	
  treeVariables3 += "nbSuperClusterBarrel/I:nbSuperClusterEndcapP/I:nbSuperClusterEndcapM/I:nbBasicClusterBarrel/I:nbBasicClusterEndcapP/I:nbBasicClusterEndcapM/I:";
  treeVariables3 +="nbSrpTTEB/I:nbZsTTEB/I:nbFlag0TTEB/I:nbFlag2TTEB/I:nbFlag4TTEB/I:nbFlag5TTEB/I:nbFlag6TTEB/I:nbFlag7TTEB/I:nbSrpTTEEM/I:nbZsTTEEM/I:nbFlag0TTEEM/I:nbFlag2TTEEM/I:nbFlag4TTEEM/I:nbFlag5TTEEM/I:nbFlag6TTEEM/I:nbFlag7TTEEM/I:nbSrpTTEEP/I:nbZsTTEEP/I:nbFlag0TTEEP/I:nbFlag2TTEEP/I:nbFlag4TTEEP/I:nbFlag5TTEEP/I:nbFlag6TTEEP/I:nbFlag7TTEEP/I:nbRH_EB/I:nbRH_EE/I";
  myEventTree_->Branch("eventTree", &(treeEvent_.eventRef),treeVariables3);

  myRecHitsTree_ = new TTree("recHitsTree","");
  TString treeVariables4 = "eventRef/I:runNum/I:bx/I:orbite/I:triggerType/I:lumiBlock/I:"; //event references
  treeVariables4 += "techTrigger0/I:techTrigger40/I:techTrigger41/I:techTrigger36/I:techTrigger37/I:techTrigger38/I:techTrigger39/I:techTrigger42/I:techTrigger43/I:"; //trigger
  treeVariables4 += "rh_barrelOrEndcap/I:rh_phi/F:rh_eta/F:rh_iPhi/I:rh_iEta/I:rh_iSm/I:rh_iSmPhi/I:rh_iSmEta/I:rh_iX/I:rh_iY/I:rh_zSide/I:rh_energy/F:rh_uncalibEnergy/F:rh_chi2/F:rh_time/F:rh_flag/I:rh_srpFlag/I:rh_isCluster/I:rh_isSuperCluster/I:"; // RecHits caracteristics
  treeVariables4 += "rh_eHfNeg/F:rh_eHfPos/F:rh_eHfNegTime/F:rh_eHfPosTime/F:rh_eHfNcounts/I:rh_eHfPcounts/I"; // HF caracteristic 
  myRecHitsTree_->Branch("recHitsTree",&(treeRecHits_.eventRef),treeVariables4);

  myTrigTree_ = new TTree("TrigTree","");
  TString treeVariables5 = "eventRef/I:runNum/I:bx/I:orbite/I:triggerType/I:lumiBlock/I:"; //event references
  treeVariables5 += "techTrigger0/I:techTrigger40/I:techTrigger41/I:techTrigger36/I:techTrigger37/I:techTrigger38/I:techTrigger39/I:"; //trigger
  treeVariables5 += "TT_barrelOrEndcap/I:TT_iphi/I:TT_ieta/I:TT_ix/I:TT_iy/I:TT_zside/I:TT_flag/I";
  myTrigTree_->Branch("TrigTree",&(treeTrigger_.eventRef),treeVariables5);

  myHFTree_ = new TTree("HFTree","");
  TString treeVariables6 = "eventRef/I:runNum/I:bx/I:orbite/I:triggerType/I:lumiBlock/I:"; //event references
  treeVariables6 += "techTrigger0/I:techTrigger40/I:techTrigger41/I:techTrigger36/I:techTrigger37/I:techTrigger38/I:techTrigger39/I:"; //trigger
  treeVariables6 += "hf_eta/F:hf_phi/F:hf_ieta/I:hf_iphi/I:hf_energy/F:hf_time/F:hf_depth/I:hf_alphaRatio/F:hf_alphaRatioTimed/F";  // HF caracteristics
  myHFTree_->Branch("HFTree",&(treeHF_.eventRef),treeVariables6);	
}
// ------------ method called once each job just after ending the event loop  ------------
void 
FirstDataAnalyzer::endJob() {
   myEventTree_->Fill(); 
   rootFile_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(FirstDataAnalyzer);
