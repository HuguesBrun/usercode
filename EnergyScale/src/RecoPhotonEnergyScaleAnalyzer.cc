// -*- C++ -*-
//
// Package:    RecoPhotonEnergyScaleAnalyzer
// Class:      RecoPhotonEnergyScaleAnalyzer
// 
/**\class RecoPhotonEnergyScaleAnalyzer RecoPhotonEnergyScaleAnalyzer.cc RecoEcal/RecoPhotonEnergyScaleAnalyzer/src/RecoPhotonEnergyScaleAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
// Original Author:  Ketino Kaadze
//         Created:  Thu Jun 21 08:59:42 CDT 2007
// $Id: RecoPhotonEnergyScaleAnalyzer.cc,v 1.2 2009/07/31 12:27:19 hbrun Exp $
//

#include "Maravin/EnergyScale/interface/RecoPhotonEnergyScaleAnalyzer.h"

//Framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TFile.h"
//Reconstruction classes
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//// Class header file
#include "RecoEcal/EgammaClusterProducers/interface/HybridClusterProducer.h"
#include "RecoEcal/EgammaCoreTools/interface/ClusterShapeAlgo.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"

#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"


#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>

#include "TString.h"
//#include "Maravin/EnergyScale/interface/fCorr_photon.h"
#include "Maravin/EnergyScale/interface/fCorr.h"
#include "Maravin/EnergyScale/interface/fEtaCorr.h"
#include "Maravin/EnergyScale/interface/AnglesUtil.h"
#include "Maravin/EnergyScale/interface/myCaloTools.h"
#include "DataFormats/Common/interface/Ref.h"

// this is from photon validator
//
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
//
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
//
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CLHEP/Units/PhysicalConstants.h"

//
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//


//========================================================================
RecoPhotonEnergyScaleAnalyzer::RecoPhotonEnergyScaleAnalyzer( const edm::ParameterSet& ps )
//========================================================================
{
  
  superClusterCollection_ = ps.getParameter<std::string>("superClusterCollection");
  superClusterProducer_   = ps.getParameter<std::string>("superClusterProducer");

  correctedSuperClusterCollection_ = ps.getParameter<std::string>("correctedSuperClusterCollection");
  correctedSuperClusterProducer_   = ps.getParameter<std::string>("correctedSuperClusterProducer");

  correctedPhotonCollection_ = ps.getParameter<std::string>("correctedPhotonCollection");

  barrelEcalHits_ = ps.getParameter<edm::InputTag>("barrelEcalHits");
  endcapEcalHits_ = ps.getParameter<edm::InputTag>("endcapEcalHits");

  hepMCLabel_ = ps.getParameter<std::string>("hepMCLabel");

  outputFile_   = ps.getParameter<std::string>("outputFile");
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms


  thePhotonMCTruthFinder_ = new PhotonMCTruthFinder();
}


//========================================================================
RecoPhotonEnergyScaleAnalyzer::~RecoPhotonEnergyScaleAnalyzer()
//========================================================================
{
  delete rootFile_;
}

//========================================================================
void
RecoPhotonEnergyScaleAnalyzer::beginJob(edm::EventSetup const&) {
//========================================================================

  mytree_ = new TTree("energyScale","");
  TString treeVariables = "mc_npar/I:parID:mc_sep/F:mc_e:mc_et:mc_phi:mc_eta:mc_theta:mc_deltaEta:mc_deltaPhi:";    // MC information
  treeVariables += "em_nEM/I:em_extra/F:em_dR/F:em_deta:em_dphi:";                 // MC <-> EM matching information
  treeVariables += "em_isInCrack/I:em_e/F:em_et:em_phi:em_eta:em_theta:em_nCell/I:em_nBC:em_r9/F:em_e5x5:em_CorrReco:"; // EM SC info    
  treeVariables	+= "em_nbcrystal/I:em_isholeinfirst/I:em_isholeinsecond/I:"; // EM SC type
  treeVariables += "pho_e/F:pho_et/F:pho_phi/F:pho_eta/F:pho_theta/F:pho_r9/F:"; // Photons variables
  treeVariables += "pho_isConverted/I:pho_nTracks/I:pho_EoverP/F:pho_Rconv/F:pho_Zconv:"; // vonvertion variable
  treeVariables += "em_pet/F:em_pe:em_peta:em_pphi:em_ptheta:"                     ;  // EM SC physics (eta corrected information)
  treeVariables += "em_phiRoad/F:em_etaRoad:"                              ;  // length of the eta/phi road of the SC

  treeVariables += "emCorrEta_e/F:emCorrEta_et/F:"                         ;  // EM f(eta) Jingzhi's correction applied
  treeVariables += "em_pw1/F:em_ew1:"                                      ;  // EM widths pw -- phiWidth, ew -- etaWidth
  treeVariables += "em_br1/F:br1_fcorr:"                                             ;  // ratios of pw/ew
  treeVariables += "emCorrBR1_e/F:emCorrBR1_et:"                           ;  // F(brem) for BR1
  treeVariables += "emCorrBR1Full_e/F:emCorrBR1Full_et:"                   ;  // Full corrections 

  treeVariables += "emCorr_e/F:emCorr_et:emCorr_eta:emCorr_phi:emCorr_theta:";// CMSSW standard corrections
  treeVariables += "emCorr_pet/F:emCorr_peta:emCorr_ptheta:"                 ;// CMSSW standard physics  

  treeVariables += "el_e/F:el_et:el_eta:el_phi:el_theta:";                    // CMSSW Electron
  treeVariables += "el_class/F:";

  treeVariables += "conv/I:conv_et/F:conv_R/F:conv_Z/F:conv_X/F:conv_Y/F:nconv/I:numb_event/I";

  mytree_->Branch("energyScale",&(tree_.mc_npar),treeVariables);


  invMassTree_ = new TTree("invMass","");
  TString massVar = "em1_e/F:em1_et:em1_eta:em1_phi:";
  massVar +=        "em2_e/F:em2_et:em2_eta:em2_phi:";
  massVar +=        "em1Corr_e/F:em1Corr_et:em1Corr_eta:em1Corr_phi:";
  massVar +=        "em2Corr_e/F:em2Corr_et:em2Corr_eta:em2Corr_phi:";
  massVar +=        "mc1_e/F:mc1_et:mc1_eta:mc1_phi:";
  massVar +=        "mc2_e/F:mc2_et:mc2_eta:mc2_phi:";
  massVar +=        "el1_e/F:el1_et:el1_eta:el1_phi:";
  massVar +=        "el2_e/F:el2_et:el2_eta:el2_phi:";
  massVar +=        "diem_mass/F:diemCorr_mass:diel_mass";                               


  // normal CMSSW di-em mass and corrected one

  invMassTree_->Branch("invMass",&(massTree_.em1_e), massVar);

}

//========================================================================
void
RecoPhotonEnergyScaleAnalyzer::analyze( const edm::Event& evt, const edm::EventSetup& es ) {
  using namespace edm; // needed for all fwk related classes
  using namespace std;
  
  //Get containers for MC truth, SC etc. ===================================================
  // =======================================================================================
  // =======================================================================================
  Handle<HepMCProduct> hepMC;
  evt.getByLabel( hepMCLabel_, hepMC ) ;
  
  const HepMC::GenEvent* genEvent = hepMC->GetEvent();
  if ( !(hepMC.isValid())) {
    LogInfo("RecoPhotonEnergyScaleAnalyzer") << "Could not get MC Product!";
    return;
  }

  //Now get the hit collection from the event ==============================================
  // =======================================================================================
  // =======================================================================================
  evt.getByLabel(barrelEcalHits_,rhcHandle_);

  if (!(rhcHandle_.isValid())) {
    edm::LogInfo("EnergyScaleAnalyzer") << "could not get a handle on the EcalRecHitCollection!";
    return;
  } 
  const EcalRecHitCollection *hit_collection = rhcHandle_.product();

  const CaloGeometry* caloGeom;
  es.get<CaloGeometryRecord>().get(theCaloGeom_);
  caloGeom = theCaloGeom_.product();

  es.get<CaloTopologyRecord>().get(theCaloTopology_);
  const CaloTopology *topology = theCaloTopology_.product();


  // get the MC truth for photons ==================================
  std::vector<SimTrack> theSimTracks;
  std::vector<SimVertex> theSimVertices;
  
  edm::Handle<SimTrackContainer> SimTk;
  edm::Handle<SimVertexContainer> SimVtx;
  evt.getByLabel("g4SimHits", SimTk);
  evt.getByLabel("g4SimHits", SimVtx);

  theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
  theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
  std::vector<PhotonMCTruth> mcPhotons=thePhotonMCTruthFinder_->find (theSimTracks,  theSimVertices);  

  //=======================Basic Cluster info
  //edm::Handle<reco::BasicClusterCollection> basicClusterHandle;
  //evt.getByLabel("hybridSuperClusters","", basicClusterHandle);
  //const reco::BasicClusterCollection* hybridBasicClusters = basicClusterHandle.product();
  
  //=======================For Vertex correction
  std::vector< Handle< HepMCProduct > > evtHandles ;
  evt.getManyByType( evtHandles ) ;
  
  for ( unsigned int i=0; i<evtHandles.size(); ++i) {
    if ( evtHandles[i].isValid() ) {
      const HepMC::GenEvent* evt = evtHandles[i]->GetEvent() ;
      
      // take only 1st vertex for now - it's been tested only of PGuns...
      //
      HepMC::GenEvent::vertex_const_iterator vtx = evt->vertices_begin() ;
      if ( evtHandles[i].provenance()->moduleLabel() == hepMCLabel_ ) {
	//Corrdinates of Vertex w.r.o. the point (0,0,0)
	xVtx_ = 0.1*(*vtx)->position().x();      
	yVtx_ = 0.1*(*vtx)->position().y();
	zVtx_ = 0.1*(*vtx)->position().z();  
      }
    }
  }
  //==============================================================================
  //Get  super clusters (normal)
  Handle<reco::SuperClusterCollection> pSuperClusters;
  try {
    evt.getByLabel(superClusterProducer_, superClusterCollection_, pSuperClusters);
  } catch ( cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer") 
      << "L210 Error! can't get collection with label " 
      << superClusterCollection_.c_str() ;
  }

  //Get CMSSW Photons
  Handle<reco::PhotonCollection> pPhotons;
  try {
    evt.getByLabel(correctedPhotonCollection_, pPhotons);
  } catch (cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer")
      << "L210 Error! can't get collection with label "
      << correctedPhotonCollection_.c_str() ;
  }

  

  // CMSSW corrected superclusters
  Handle<reco::SuperClusterCollection> pCorrectedSuperClusters;
  try {
    evt.getByLabel(correctedSuperClusterProducer_, 
		   correctedSuperClusterCollection_, 
		   pCorrectedSuperClusters);
  } catch ( cms::Exception& ex ) {
    edm::LogError("RecoPhotonEnergyScaleAnalyzer") 
      << "L221 Error! can't get collection with label " 
      << correctedSuperClusterCollection_.c_str() ;
  }

  
  
  
  // All containers are loaded, perform the analysis///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin();
  
  // Search for MC electrons or photons that satisfy the criteria
  float min_eT = 1.9;
  float max_eta = 1.6;
  
  std::vector<HepMC::GenParticle* > mcParticles;
  //int counter = 0;
  for ( HepMC::GenEvent::particle_const_iterator p = genEvent->particles_begin(); 
	p != genEvent->particles_end(); 
	++p ) {
    //LogInfo("RecoPhotonEnergyScaleAnalyzer") << "Particle " << ++counter 
    //<< " PDG ID = " << (*p)->pdg_id() << " pT = " << (*p)->momentum().perp();
    // require photon or electron
    if ( (*p)->pdg_id() != 22 && abs((*p)->pdg_id()) != 11 ) continue;
    
    // require selection criteria
    bool satisfySelectionCriteria = 
      (*p)->momentum().perp() > min_eT &&
      fabs((*p)->momentum().eta()) < max_eta;
    
    if ( ! satisfySelectionCriteria ) continue;

    // EM MC particle is found, save it in the vector
    mcParticles.push_back(*p);
  }
  // separation in dR between 2 first MC particles
  // should not be used for MC samples with > 2 em objects generated!
  if ( mcParticles.size() == 2 ) {
    HepMC::GenParticle* mc1 = mcParticles[0];
    HepMC::GenParticle* mc2 = mcParticles[1];
    tree_.mc_sep = kinem::delta_R(mc1->momentum().eta(), mc1->momentum().phi(),
				  mc2->momentum().eta(), mc2->momentum().phi());
  } else
    tree_.mc_sep = -100;

  /*  // Quick fix for double pair of electrons ====================================
  vector<float> etas;
  vector<float> phis;
  for(std::vector<HepMC::GenParticle* >::const_iterator p = mcParticles.begin();
      p != mcParticles.end(); ++p) {
    HepMC::GenParticle* mc = *p;
    std::vector<float>::iterator ie = etas.end();
    float eta = fabs(mc->momentum().eta());
    float phi = fabs(mc->momentum().phi());
    ie = find(etas.begin(), etas.end(), eta);
    if ( ie == etas.end() ) {
      etas.push_back(eta);
      phis.push_back(phi);
    }
  }
  
  if (etas.size() == 2 ) {
    // found 4 electrons!
    tree_.mc_deltaEta = fabs(etas[0]-etas[1]);
    tree_.mc_deltaPhi = kinem::delta_phi(phis[0], phis[1]);
    if ( tree_.mc_deltaPhi > kinem::PI/2 ) 
      tree_.mc_deltaPhi = kinem::PI - tree_.mc_deltaPhi;
  }
  else                   {
    tree_.mc_deltaEta = 100.;
    tree_.mc_deltaPhi = 100.;
    }*/
  // remove all the events with a potential to mix 4 leptons in the final state!
  //if ( tree_.mc_deltaEta < 0.4 ) continue;
  
  // now loop over MC particles, find the match with SC and do everything we need
  // then save info in the tree for every MC particle
  std::vector<math::XYZTLorentzVector> p4;
  std::vector<math::XYZTLorentzVector> p4Corr;
  std::vector<math::XYZTLorentzVector> p4mc;
  std::vector<math::XYZTLorentzVector> p4el;

  int nEM = 0;
  for(std::vector<HepMC::GenParticle* >::const_iterator p = mcParticles.begin();
      p != mcParticles.end(); ++p) {
    nEM = 0;
    HepMC::GenParticle* mc = *p;
    
    // Fill MC information
    tree_.mc_npar  = mcParticles.size();
    tree_.numb_event = evt.id().event();	
    tree_.parID    = mc->pdg_id();
    tree_.mc_e     = mc->momentum().e();
    tree_.mc_et    = mc->momentum().e()*sin(mc->momentum().theta());
    tree_.mc_phi   = mc->momentum().phi();
    tree_.mc_eta   = mc->momentum().eta();
    tree_.mc_theta = mc->momentum().theta();
    
    const reco::SuperClusterCollection* superClusters = pSuperClusters.product();
    
    //Loop over super clusters and pick up the most energetic supercluster
    //that is matched with MC particle mc (dR < 0.2)
    
    reco::SuperClusterCollection::const_iterator em = superClusters->end();
    float energyMax = -100.0; // dummy energy of the matched SC
    int nbSChere = 0;
    for(reco::SuperClusterCollection::const_iterator aClus = superClusters->begin();
	aClus != superClusters->end(); ++aClus) {
      // check the matching
      nbSChere++;	
      float dR = kinem::delta_R(mc   ->momentum().eta(), mc   ->momentum().phi(), 
				aClus->position().eta(), aClus->position().phi());
  //    std::cout << "le dR = " << dR << std::endl;
      if (dR <  1.0) { // a rather loose matching cut
	float energy = aClus->energy();
//	std::cout << "on a energy = " << energy << std::endl;
	if ( energy < energyMax ) continue;
//	std::cout << "on garde energy = " << energy << std::endl;
	energyMax = energy;
	em = aClus;
	++nEM; // how many superclusters are associated with this electron?
      }
    }
    if (nbSChere > 2 ) std::cout << " Ya plus de deux SC !!!! " << std::endl;
    tree_.em_nEM = nEM;

    if (em == superClusters->end() ) {
	continue; // no matching EM object is found, go to the next MC particle
	std::cout << "pas de maching avec un super cluster " << std::endl;
    }
//    std::cout << "ya bien un matching !!! " << std::endl;

    tree_.em_extra = 0;
    for(reco::SuperClusterCollection::const_iterator aClus = superClusters->begin();
        aClus != superClusters->end(); ++aClus) {
      if ( em == aClus ) continue;

      // check the matching
      float dR = kinem::delta_R(mc   ->momentum().eta(), mc   ->momentum().phi(),
                                aClus->position().eta(), aClus->position().phi());
      if (dR <  1.0) { // a rather loose matching cut
        tree_.em_extra += aClus->energy();
      }
    }


    //looking for the photon associated with the supercluster
    reco::PhotonCollection::const_iterator myphoton;  
    const reco::PhotonCollection* photons = pPhotons.product();
    int nbphoton = 0;
    for(reco::PhotonCollection::const_iterator aPho = photons->begin(); aPho != photons->end(); aPho++){
      const reco::SuperCluster *theSC = aPho->superCluster().get();
      if (theSC->position().eta() == em->position().eta()){
	myphoton = aPho;
	nbphoton++;
      }
    }
	
    if (nbphoton>0){
	tree_.pho_e = myphoton->energy();
	tree_.pho_et = myphoton->energy() * sin(myphoton->theta());
	tree_.pho_r9 = myphoton->r9();
	tree_.pho_eta = myphoton->eta();
	tree_.pho_phi = myphoton->phi();
	tree_.pho_theta = myphoton->theta();
	reco::ConversionRefVector convrefvector = myphoton->conversions();
	if (convrefvector.size() == 1){
	    reco::ConversionRef lareference = convrefvector[0];
	    const reco::Conversion *laconversion = lareference.get();
	    tree_.pho_isConverted = 1;
	    tree_.pho_EoverP = laconversion->EoverP();
	    tree_.pho_Rconv = sqrt(laconversion->conversionVertex().x()*laconversion->conversionVertex().x()+laconversion->conversionVertex().y()*laconversion->conversionVertex().y());
            tree_.pho_Zconv = laconversion->conversionVertex().z();   
	}
	else {
            tree_.pho_isConverted = 0;
            tree_.pho_EoverP = -99;
            tree_.pho_Rconv = -99;
            tree_.pho_Zconv = -99;
	}
    }
    else {
//	std::cout << "Pas de photon associé à ce Super Cluster !!! " << std::endl;
	   tree_.pho_e  = -100; 
    }	
    // Matching SC found, add info to the tree
    tree_.em_isInCrack = 0;
    double emAbsEta = fabs(em->position().eta());
    // copied from RecoEgama/EgammaElectronAlgos/src/EgammaElectronClassification.cc
    if ( emAbsEta < 0.018 ||
	 (emAbsEta > 0.423 && emAbsEta < 0.461) || 
	 (emAbsEta > 0.770 && emAbsEta < 0.806) || 
	 (emAbsEta > 1.127 && emAbsEta < 1.163) || 
	 (emAbsEta > 1.460 && emAbsEta < 1.558) )
      tree_.em_isInCrack = 1;
    
    tree_.em_dR = kinem::delta_R(mc->momentum().eta(), mc->momentum().phi(), 
				 em->position().eta(), em->position().phi());
    tree_.em_deta = (mc->momentum().eta() - em->position().eta());
    tree_.em_dphi = kinem::delta_phi(mc->momentum().phi(), em->position().phi());
    tree_.em_e5x5 = EcalClusterTools::e5x5( *(em->seed() ), hit_collection, &(*topology));
    tree_.em_e     = em->rawEnergy();
//    std::cout << " energie RAW du SC = " << em->rawEnergy() << std::endl;
    tree_.em_et    = em->rawEnergy() * sin(mc->momentum().theta());
    tree_.em_phi   = em->position().phi();
    tree_.em_eta   = em->position().eta();
    tree_.em_theta = em->position().theta();
    tree_.em_CorrReco = em->energy();
    tree_.em_nCell = 0;//(em->getHitsByDetId()).size();
    tree_.em_nBC   = em->clustersSize();
    if (em->clustersSize() == 1 ) {
	const reco::BasicCluster leseed = *(em->seed());
	const std::vector< std::pair<DetId, float> > listcristaux = leseed.hitsAndFractions();
	tree_.em_nbcrystal = listcristaux.size();
	DetId seedXtalId = EcalClusterTools::getMaximum( leseed.hitsAndFractions(), hit_collection).first;
	//std::cout << myCaloTools::deltaIphi(359,1) << std::endl;
	tree_.em_isholeinfirst = myCaloTools::contourinf(leseed,seedXtalId);
	tree_.em_isholeinsecond = myCaloTools::contoursup(leseed,seedXtalId);
//	std::cout << "nbXtal = " << listcristaux.size() << " nb petite couronne= " << bin1 << " nb grande couronne = " << myCaloTools::contoursup(leseed,seedXtalId) << std::endl;

    } 
    else{
	tree_.em_nbcrystal = -10;
	tree_.em_isholeinfirst = -10;
	tree_.em_isholeinsecond = -10;
    }

    float e3x3 = EcalClusterTools::e3x3( *(em->seed() ), hit_collection, &(*topology));
    float r9 = e3x3/em->rawEnergy();
    tree_.em_r9 = r9;   


    
    //Get physics e, et etc:
    //Coordinates of EM object with respect of the point (0,0,0)
    xClust_zero_ = em->position().x();
    yClust_zero_ = em->position().y();
    zClust_zero_ = em->position().z();
    //Coordinates of EM object w.r.o. the Vertex position
    xClust_vtx_ = xClust_zero_ - xVtx_;
    yClust_vtx_ = yClust_zero_ - yVtx_;
    zClust_vtx_ = zClust_zero_ - zVtx_;
    
    energyMax_ = em->energy();
    thetaMax_ = em->position().theta();
    etaMax_ = em->position().eta();
    phiMax_ = em->position().phi();
    eTMax_ = energyMax_ * sin(thetaMax_);
    if ( phiMax_ < 0) phiMax_ += 2*3.14159;
    
    rClust_vtx_ = sqrt(xClust_vtx_*xClust_vtx_ + yClust_vtx_*yClust_vtx_ + zClust_vtx_*zClust_vtx_);
    thetaMaxVtx_ = acos(zClust_vtx_/rClust_vtx_);
    etaMaxVtx_   = - log(tan(thetaMaxVtx_/2));
    eTMaxVtx_    = energyMax_ * sin(thetaMaxVtx_); 
    phiMaxVtx_   = atan2(yClust_vtx_,xClust_vtx_); 
    if ( phiMaxVtx_ < 0) phiMaxVtx_ += 2*3.14159;
    //=============================
    //parametres of EM object after vertex correction
    tree_.em_pet    = eTMaxVtx_;
    tree_.em_pe     = tree_.em_pet/sin(thetaMaxVtx_);
  /*  tree_.em_peta   = myphoton->eta();
    tree_.em_pphi = myphoton->phi();	*/
    tree_.em_ptheta = thetaMaxVtx_;

    //========================================================================================
    // Information for the dynamic road parameters
    float phiMin = 100.;
    float phiMax = -100.;
    float etaMin = 100.;
    float etaMax = -100.;
    /*
    //get recHit ID which comprise super cluster
    std::vector<DetId> detId = em->getHitsByDetId();
    
    for (std::vector<DetId>::iterator hit = detId.begin();
    hit != detId.end(); ++hit) { 
    // Loop over recHits associated with the given supercluster
      
    EcalRecHitCollection::const_iterator rHit = hit_collection->find(*hit);
    // rHit will point to the recHit in the ECAL recHit collection
    //Get the cell geometry of a given detector id 
    const CaloCellGeometry *this_cell = (geometry).getGeometry(rHit->id());  
    if ( this_cell == 0 ) {
    edm::LogInfo("RecoPhotonEnergyScaleAnalyzer") << "pointer to the cell is NULL!"; 
    return;
    }
    GlobalPoint position = this_cell->getPosition();  //returns the position of reference for this cell
    float phiHit = position.phi();
    float etaHit = position.eta();
    if ( phiHit < phiMin ) phiMin = phiHit;
    if ( phiHit > phiMax ) phiMax = phiHit;
    if ( etaHit < etaMin ) etaMin = etaHit;
    if ( etaHit > etaMax ) etaMax = etaHit;
    }
    */

    tree_.em_phiRoad = kinem::delta_phi(phiMin, phiMax);
    tree_.em_etaRoad = kinem::delta_eta(etaMax, etaMin);
    //cout << tree_.em_phiRoad << " " << tree_.em_etaRoad << endl;

    float absEta = fabs(em->position().eta());
    float emTheta = em->position().theta();

    float energyWithEtaCorrection = tree_.em_e/fcorr::f5x5((absEta*(5/0.087)));
    
    float pTWithEtaCorrection     = energyWithEtaCorrection * sin(emTheta);
    tree_.emCorrEta_e = energyWithEtaCorrection;
    tree_.emCorrEta_et = pTWithEtaCorrection;
    
    tree_.em_pw1 = em->phiWidth();
    tree_.em_ew1 = em->etaWidth();
    //    std::cout << " phi = " << em->phiWidth() << " eta = " << em->etaWidth() << std::endl; 
    tree_.em_br1 = tree_.em_pw1/tree_.em_ew1;
    
    // Corrected energies
    tree_.emCorrBR1_e  = fcorr::electron_br1(tree_.em_br1, tree_.emCorrEta_e);
	
    tree_.br1_fcorr = tree_.emCorrEta_e/tree_.emCorrBR1_e;
    tree_.emCorrBR1_et = tree_.emCorrBR1_e * sin(em->position().theta());
    
    tree_.emCorrBR1Full_et = fcorr::electron_br1_complete(tree_.emCorrBR1_et, absEta);
    tree_.emCorrBR1Full_e  = tree_.emCorrBR1Full_et/sin(emTheta);


    // =======================================================================
    // Get conversion parameters
    // =======================================================================

    tree_.conv = 0;
    tree_.nconv = 0;

    for ( std::vector<PhotonMCTruth>::const_iterator mcPho = mcPhotons.begin(); mcPho != mcPhotons.end(); mcPho++) {
      // is it matched?
      double eta = (*mcPho).fourMomentum().eta();
      double phi = (*mcPho).fourMomentum().phi();
      
      double dR = kinem::delta_R(eta, phi, mc->momentum().eta(), mc->momentum().phi());
      if ( dR < 1e-6 ) {
	tree_.conv = (*mcPho).isAConversion();
	tree_.conv_et = (*mcPho).fourMomentum().et();
	tree_.conv_R = (*mcPho).vertex().perp();
	tree_.conv_Z = (*mcPho).vertex().mag() * (*mcPho).vertex().cosTheta();
	tree_.conv_X = (*mcPho).vertex().x();
	tree_.conv_Y = (*mcPho).vertex().y();
      }
      if ( dR < 0.1 ) tree_.nconv += 1;
    } 
    
    /*
    // =======================================================================
    // Get PixelMatchGsfElectrons
    // =======================================================================
    const reco::PixelMatchGsfElectronCollection* electrons = pGsfElectrons.product();
    for(reco::PixelMatchGsfElectronCollection::const_iterator iel = electrons->begin();
    iel != electrons->end(); ++iel ) {
    // check the matching
    double dR = kinem::delta_R(mc->momentum().eta(), mc->momentum().phi(),
    iel->p4().eta(), iel->p4().phi());
    if ( dR < 0.7 ) {
    tree_.el_e  = iel->p4().e();
    tree_.el_et = iel->p4().Et();
    tree_.el_eta = iel->p4().eta();
    tree_.el_phi = iel->p4().phi();
    tree_.el_theta = iel->p4().theta();
    tree_.el_class = iel->classification();
    math::XYZTLorentzVector p4(tree_.el_et*cos(tree_.el_phi), 
    tree_.el_et*sin(tree_.el_phi), 
    tree_.el_e*cos(tree_.el_theta), 
    tree_.el_e);
    //std::cout << iel->p4().Et() << " " << iel->gsfTrack()->outerMomentum().Rho() << std::endl;
    p4el.push_back(p4);
    }
    }
    */
    
    // =======================================================================
    // Get super clusters after energy correction (Standard CMSSW corrections)
    // =======================================================================
    
    const reco::SuperClusterCollection* correctedSuperClusters = pCorrectedSuperClusters.product();
    //Loop over super clusters and pick up the most energetic supercluster
    //that is matched with MC particle mc (dR < 0.2)
    
    em = correctedSuperClusters->end();
    energyMax = -100.0; // dummy energy of the matched SC
    for(reco::SuperClusterCollection::const_iterator aClus = correctedSuperClusters->begin();
	aClus != correctedSuperClusters->end(); ++aClus) {
      // check the matching
      float dR = kinem::delta_R(mc   ->momentum().eta(), mc   ->momentum().phi(), 
				aClus->position().eta(), aClus->position().phi());
      if (dR <  0.7) {
	float energy = aClus->energy();
	if ( energy < energyMax ) continue;
	energyMax = energy;
	em = aClus;
      }
    }
    if (em == correctedSuperClusters->end() ) continue; // no matching EM object is found, go to the next MC particle
    
    //fill tree with kinematic variables of corrected Super Cluster
    tree_.emCorr_e     = em->energy();
    tree_.emCorr_et    = em->energy() * sin(em->position().theta());
    tree_.emCorr_phi   = em->position().phi();
    tree_.emCorr_eta   = em->position().eta();
    tree_.emCorr_theta = em->position().theta();
    
    //===========================Corrects vertex position for most energetic cluster after energy correction is applied
    //Coordinates of EM object with respect of the point (0,0,0)
    xClust_zero_ = em->position().x();
    yClust_zero_ = em->position().y();
    zClust_zero_ = em->position().z();
    //Coordinates of EM object w.r.o. the Vertex position
    xClust_vtx_ = xClust_zero_ - xVtx_; 
    yClust_vtx_ = yClust_zero_ - yVtx_;
    zClust_vtx_ = zClust_zero_ - zVtx_;
    
    energyMax_ = em->energy();  
    thetaMax_ = em->position().theta();
    etaMax_ = em->position().eta();
    phiMax_ = em->position().phi();
    eTMax_ = energyMax_ * sin(thetaMax_);
    if ( phiMax_ < 0) phiMax_ += 2*3.14159;
    
    rClust_vtx_ = sqrt (xClust_vtx_*xClust_vtx_ + yClust_vtx_*yClust_vtx_ + zClust_vtx_*zClust_vtx_);
    thetaMaxVtx_ = acos(zClust_vtx_/rClust_vtx_);
    etaMaxVtx_   = - log(tan(thetaMaxVtx_/2));
    eTMaxVtx_    = energyMax_ * sin(thetaMaxVtx_); 
    phiMaxVtx_   = atan2(yClust_vtx_,xClust_vtx_); 
    if ( phiMaxVtx_ < 0) phiMaxVtx_ += 2*3.14159;
    //================
    //parametres of EM object after vertex correction (energy correction is applied)
    tree_.emCorr_pet    = eTMaxVtx_;
    tree_.emCorr_peta   = etaMaxVtx_;
    tree_.emCorr_ptheta = thetaMaxVtx_;
    
    // Filling the LorentzVectors ====================================================
    // This is corrected em Lorentz vector
    double e = tree_.emCorrBR1Full_e;
    double px = tree_.emCorrBR1Full_et*cos(tree_.em_phi);
    double py = tree_.emCorrBR1Full_et*sin(tree_.em_phi);
    double pz = tree_.emCorrBR1Full_e*cos(tree_.em_ptheta);
    math::XYZTLorentzVector emCorr_p4(px, py, pz, e);
    p4Corr.push_back(emCorr_p4);
    // This is standard CMSSW
    e = tree_.emCorr_pet*cosh(tree_.emCorr_peta);
    px = tree_.emCorr_pet*cos(tree_.em_phi);
    py = tree_.emCorr_pet*sin(tree_.em_phi);
    pz = e*cos(tree_.em_ptheta);
    math::XYZTLorentzVector em_p4(px, py, pz, e);
    p4.push_back(em_p4);
    // This is generate MC
    e = tree_.mc_e;
    px = tree_.mc_et*cos(tree_.mc_phi);
    py = tree_.mc_et*sin(tree_.mc_phi);
    pz = e*cos(tree_.mc_theta);
    math::XYZTLorentzVector mc_p4(px, py, pz, e);
    p4mc.push_back(mc_p4);
 
    mytree_->Fill();
  } // loop over particles
  
  
    // Produce invariant mass variables
  if ( p4Corr.size() > 1 && p4.size() > 1 && p4el.size() > 1 ) {
    math::XYZTLorentzVector diemCorr(p4Corr[0] + p4Corr[1]);
    math::XYZTLorentzVector diem    (p4[0]     + p4[1]);
    math::XYZTLorentzVector diel    (p4el[0]   + p4el[1]);
    
    massTree_.em1_e     = p4[0].e();
    massTree_.em1_et    = p4[0].Pt();
    massTree_.em1_eta   = p4[0].eta();
    massTree_.em1_phi   = p4[0].phi();
    massTree_.em2_e     = p4[1].e();
    massTree_.em2_et    = p4[1].Pt();
    massTree_.em2_eta   = p4[1].eta();
    massTree_.em2_phi   = p4[1].phi();
    
    massTree_.em1Corr_e     = p4Corr[0].e();
    massTree_.em1Corr_et    = p4Corr[0].Pt();
    massTree_.em1Corr_eta   = p4Corr[0].eta();
    massTree_.em1Corr_phi   = p4Corr[0].phi();
    massTree_.em2Corr_e     = p4Corr[1].e();
    massTree_.em2Corr_et    = p4Corr[1].Pt();
    massTree_.em2Corr_eta   = p4Corr[1].eta();
    massTree_.em2Corr_phi   = p4Corr[1].phi();
    
    massTree_.el1_e   = p4el[0].e();
    massTree_.el1_et  = p4el[0].Pt();
    massTree_.el1_eta = p4el[0].eta();
    massTree_.el1_phi = p4el[0].phi();
    massTree_.el1_e   = p4el[1].e();
    massTree_.el1_et  = p4el[1].Pt();
    massTree_.el1_eta = p4el[1].eta();
    massTree_.el1_phi = p4el[1].phi();
    
    massTree_.diem_mass     = diemCorr.M2();
    massTree_.diemCorr_mass = diem.    M2();
    massTree_.diel_mass     = diel.    M2();
    
    massTree_.mc1_e   = p4mc[0].e();
    massTree_.mc1_et  = p4mc[0].Pt();
    massTree_.mc1_eta = p4mc[0].eta();
    massTree_.mc1_phi = p4mc[0].phi();
    
    massTree_.mc2_e   = p4mc[1].e();
    massTree_.mc2_et  = p4mc[1].Pt(); 
    massTree_.mc2_eta = p4mc[1].eta();
    massTree_.mc2_phi = p4mc[1].phi();
    
    invMassTree_->Fill();
  }
  
}

//==========================================================================
void
RecoPhotonEnergyScaleAnalyzer::endJob() {
  //========================================================================
  //Fill ROOT tree
  rootFile_->Write();
}

DEFINE_FWK_MODULE(RecoPhotonEnergyScaleAnalyzer);
