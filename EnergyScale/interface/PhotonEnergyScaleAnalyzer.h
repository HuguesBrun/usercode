#ifndef Maravin_EgammaClusterProducers_PhotonEnergyScaleAnalyzer_h
#define Maravin_EgammaClusterProducers_PhotonEnergyScaleAnalyzer_h
/**\class PhotonEnergyScaleAnalyzerx

 Description: Analyzer to fetch collection of objects from event and make simple plots

 Implementation:
     \\\author: Keti Kaadze, June 2007
*/
//
// $Id: PhotonEnergyScaleAnalyzer.h,v 1.1.1.1 2009/07/10 14:09:15 hbrun Exp $
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//Geometry
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

// Ecal rec hits
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TNtuple.h"
class TFile;

//
// class declaration
//

class PhotonEnergyScaleAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PhotonEnergyScaleAnalyzer( const edm::ParameterSet& );
      ~PhotonEnergyScaleAnalyzer();


      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void beginJob(edm::EventSetup const&);
      virtual void endJob();
 private:

      edm::ESHandle<CaloGeometry> theCaloGeom_;
      edm::ESHandle<CaloTopology> theCaloTopology_;

      edm::Handle<EcalRecHitCollection> rhcHandle_;

      edm::InputTag barrelEcalHits_;
      edm::InputTag endcapEcalHits_;

      PhotonMCTruthFinder*  thePhotonMCTruthFinder_;


      std::string outputFile_; // output file
  
      std::string superClusterCollection_;
      std::string superClusterProducer_;
      
      std::string correctedSuperClusterCollection_;
      std::string correctedSuperClusterProducer_;
      
      std::string hepMCLabel_;

      std::string hitProducer_;
      std::string hitCollection_;
      
      // root file to store histograms
      TFile*  rootFile_;

      // InvMass tree
      TTree* invMassTree_;

      struct invMassTree_structure_ {
	float em1_e;
	float em1_et;
	float em1_eta;
	float em1_phi;
	float em2_e;
	float em2_et;
	float em2_eta;
	float em2_phi;

	float em1Corr_e;
	float em1Corr_et;
	float em1Corr_eta;
	float em1Corr_phi;
	float em2Corr_e;
	float em2Corr_et;
	float em2Corr_eta;
	float em2Corr_phi;

	float el1_e;
	float el1_et;
	float el1_eta;
	float el1_phi;
	float el2_e;
	float el2_et;
	float el2_eta;
	float el2_phi;

	float mc1_e;
	float mc1_et;
	float mc1_eta;
	float mc1_phi;
	float mc2_e;
	float mc2_et;
	float mc2_eta;
	float mc2_phi;

	float diem_mass;
	float diemCorr_mass;
	float diel_mass;
      };
      invMassTree_structure_ massTree_;

      //Tree
      TTree* mytree_;
      struct tree_structure_ {
	// MC information
	int   mc_npar;
	int   parID;
	float mc_sep;
	float mc_e;
	float mc_et;
	float mc_phi;
	float mc_eta;
	float mc_theta;
	float mc_deltaEta;
	float mc_deltaPhi;

	// MC-EM matching info
	int   em_nEM;
	float em_extra;
	float em_dR;
	float em_deta;
	float em_dphi;
	
	// EM SC info (uncorrected)
	int   em_isInCrack;
	float em_e;
	float em_et;
	float em_phi;
	float em_eta;
	float em_theta;
	int   em_nCell;
	int   em_nBC;
	float em_r9;
	float em_e5x5;

	// physics variables
	float em_pet;
	float em_pe;
	float em_peta;
	float em_ptheta;

	float em_phiRoad; // number of crystals of the SC in phi
	float em_etaRoad; // -//- the same in eta

	// EM f(eta) Jingzhi's correction
	float emCorrEta_e;
	float emCorrEta_et;
	
	// EM widths, pw -- phiWidth, ew -- etaWidth
	// 1 -- linear, 2 -- log(e/0.1), 3 -- log(E/Esc + 1) weighting
	float em_pw1;
	float em_ew1;

	// ratios of widths pw/ew
	float em_br1;
	float br1_fcorr;

	// F(brem) corrections
	float emCorrBR1_e;
	float emCorrBR1_et;

	// Full corrections
	float emCorrBR1Full_e;
	float emCorrBR1Full_et;

	// CMSSW corrections
	float emCorr_e;
	float emCorr_et;
	float emCorr_eta;
	float emCorr_phi;
	float emCorr_theta;
	float emCorr_pet;
	float emCorr_peta;
	float emCorr_ptheta;

	// CMSSW electrons
	float el_e;
	float el_et;
	float el_eta;
	float el_phi;
	float el_theta;
	float el_class;

	// CMSSW photon conversion (MC Truth info)
	int conv;
	float conv_et;
	float conv_R;
	int nconv;
      };
      tree_structure_ tree_;
      //
      float xVtx_;
      float yVtx_;
      float zVtx_;
      //
      float xClust_zero_;
      float yClust_zero_;
      float zClust_zero_;
      //
      float xClust_vtx_;
      float yClust_vtx_;
      float zClust_vtx_;
      //
      float rClust_vtx_;
      //
      float energyMax_;
      float eTMax_;
      float eTMaxVtx_;
      float etaMax_;
      float etaMaxVtx_;
      float phiMax_;
      float phiMaxVtx_;
      float thetaMax_;
      float thetaMaxVtx_;
      //
};
#endif

