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
// $Id: FirstDataAnalyzer.h,v 1.3 2010/04/28 09:12:52 hbrun Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>


#include "TString.h"
#include "TFile.h"
#include "TTree.h" 

// user include files
#include "myCaloTools.h"
#include "fEtaCorr.h"
#include "fCorr.h"
#include "fCorr_EE.h"
#include "AnglesUtil.h"

// FrameWork include files ..........................................................

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

// CMS detector geometry and topology and calibration 
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"


// L1 trigger 
#include "DataFormats/L1GlobalTrigger/interface/L1GtTechnicalTrigger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// HLT trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

// HF rechits
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

// REchit collection
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h" 

// SRP RecHits flags collection
#include "DataFormats/EcalDigi/interface/EBSrFlag.h"
#include "DataFormats/EcalDigi/interface/EESrFlag.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Trigger Tower Map collection and DetId
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "DataFormats/EcalDetId/interface/EcalScDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"

// Basic Cluster collection
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

// Super Cluster Collection
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// Photons Collection
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
// Converted photon
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
// Electrons Collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//MC truth
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"

// The SC Tools
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"




// class decleration
//
//class TFile;

class FirstDataAnalyzer : public edm::EDAnalyzer {
   public:
      explicit FirstDataAnalyzer(const edm::ParameterSet&);
      ~FirstDataAnalyzer();


   private:
      virtual void mySuperClusterAnalyzer(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::InputTag triggerL1Tag_, edm::InputTag HFrecHitsCollection_, edm::InputTag ecalHits_, std::string srpFlagsCollection_, std::string clusterCollection_,std::string clusterProducer_,std::string correctedSuperClusterCollection_,std::string correctedSuperClusterProducer_,std::string photonCollection_, std::string electronCollection_, std::string MCParticlesCollection_, bool isBarrel);
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


// Global Handles
//              Geometry and Topology of the detector
      edm::ESHandle<CaloGeometry> theCaloGeom_;
      edm::ESHandle<CaloTopology> theCaloTopology_;




//              Electromagnetic Calo Rec hits collections
/*      edm::Handle<EcalRecHitCollection> rhcHandleBarrel_;
      edm::Handle<EcalRecHitCollection> rhcHandleEndcap_;*/

      int nbSCEB_;
      int nbSCEEP_;
      int nbSCEEM_;
      int nbBCEB_;
      int nbBCEEP_;
      int nbBCEEM_;

  // Selective ReadOut State
      int nbSrpTTEB_;
      int nbZsTTEB_;
      int nbFlag0TTEB_;
      int nbFlag2TTEB_;
      int nbFlag4TTEB_;
      int nbFlag5TTEB_;
      int nbFlag6TTEB_;
      int nbFlag7TTEB_;
      int nbSrpTTEEM_;
      int nbZsTTEEM_;
      int nbFlag0TTEEM_;
      int nbFlag2TTEEM_;
      int nbFlag4TTEEM_;
      int nbFlag5TTEEM_;
      int nbFlag6TTEEM_;
      int nbFlag7TTEEM_;
      int nbSrpTTEEP_;
      int nbZsTTEEP_;
      int nbFlag0TTEEP_;
      int nbFlag2TTEEP_;
      int nbFlag4TTEEP_;
      int nbFlag5TTEEP_;
      int nbFlag6TTEEP_;
      int nbFlag7TTEEP_;
      int nbRH_EB_;
      int nbRH_EE_;
	


// Read the python config files

//              Do I put the Rec Hits collection in the Tree ?
      bool isDoRecHits_;
//		Do I put the TT infos ?
      bool isDoTTflag_;	
//		Do I put the HF recHits ? 
      bool isDoHFrecHits_; 	
//  		Do I read MC truth ?
      bool isMCTruth_;
//		Do I read photons MC truth ?
      bool isPhotonMCTruth_;
//		Delta R for Matching with MC truth !
      double deltaRMax_;	
//		Trigger L1 technical bytes
      edm::InputTag triggerL1Tag_;
//		HTL trigger byte
      edm::InputTag triggerHLTTag_;		
// 		Forward Hadronique Calorimeter recHits
      edm::InputTag HFrecHitsCollection_;
//              Electromagnetic Calo Rec hits collections
      edm::InputTag barrelEcalHits_;
      edm::InputTag endcapEcalHits_;
//		ECAL SRP Flags collection
      std::string barrelSrpFlagsCollection_;
      std::string endcapSrpFlagsCollection_;
//		ECAL ClusterCollections
      std::string barrelClusterCollection_;
      std::string barrelClusterProducer_;
      std::string endcapClusterCollection_;
      std::string endcapClusterProducer_;
//		ECAL Super Cluster Collections
      std::string barrelCorrectedSuperClusterCollection_;
      std::string barrelCorrectedSuperClusterProducer_;
      std::string endcapCorrectedSuperClusterCollection_;
      std::string endcapCorrectedSuperClusterProducer_;
//		Photons Collections
      std::string photonCollection_;
//		Electrons Collections
      std::string electronCollection_;	
//		MC particles collections
      std::string MCParticlesCollection_;
//		Geant 4 sim hits 
      std::string Geant4SimHitsCollection_;
//  		Name of the output File
      std::string outputFile_;
	

// the Nancy MC truth analyzer
      PhotonMCTruthFinder*  thePhotonMCTruthFinder_;


// OutFile with the output Tree
      TFile* rootFile_;

// Root Tree declaration
      TTree* myTree_;
      TTree* myBCTree_;
      TTree* myEventTree_;
      TTree* myRecHitsTree_;	
      TTree* myTrigTree_;
      TTree* myHFTree_;	    
// Declaration of all the TTree leaves in a struct
      struct tree_structure_ {
		//references to the Event
	        int eventRef;
	        int runNum;
	        int bx;
        	int orbite;
	        int triggerType;
		int lumiBlock;
	        // technical L1 trigger
	        int techTrigger0;
	        int techTrigger40;
	        int techTrigger41;
	        int techTrigger36;
	        int techTrigger37;
	        int techTrigger38;
	        int techTrigger39;
		// HTL trigger
		int HLT_Photon10_L1R;
		int HLT_Photon15_L1R;
		int HLT_Photon15_LooseEcalIso_L1R;
		int HLT_Photon20_L1R;
		int HLT_Photon30_L1R_8E29;
		int HLT_DoublePhoton4_Jpsi_L1R;
		int HLT_DoublePhoton4_Upsilon_L1R;
		int HLT_DoublePhoton4_eeRes_L1R;
		int HLT_DoublePhoton5_eeRes_L1R;
		int HLT_DoublePhoton5_Jpsi_L1R;
		int HLT_DoublePhoton5_Upsilon_L1R;
		int HLT_DoublePhoton5_L1R;
		int HLT_DoublePhoton10_L1R;
		int HLT_DoubleEle5_SW_L1R;
		int HLT_Ele20_LW_L1R;
		int HLT_Ele15_SiStrip_L1R;
		int HLT_Ele15_SC10_LW_L1R;
		int HLT_Ele15_LW_L1R;
		int HLT_Ele10_LW_EleId_L1R;
		int HLT_Ele10_LW_L1R;
		int HLT_Photon15_TrackIso_L1R;

		// info about the Super-Cluster
	        int   em_isInCrack;
		int   em_barrelOrEndcap;
	        float em_e;
		float em_eRAW;
	        float em_et;
		float em_etRAW;
	        float em_phi;
	        float em_eta;
	        float em_theta;
	        float em_e5x5;
	        float em_e2x2; 	
		// shape of the super cluster
	        float em_pw1;
   	        float em_ew1;
		float em_sigmaetaeta;
		float em_sigmaietaieta;
		float em_sigmaphiphi;
		float em_sigmaiphiiphi;
                float em_BCsigmaetaeta;
                float em_BCsigmaietaieta;
                float em_BCsigmaphiphi;
                float em_BCsigmaiphiiphi;
                float em_r9;
		float em_r19;
	        float em_br1;
                int   em_nBC;
		int   em_nbcrystal;
	        int   em_nbcrystalSEED;
	        int   em_isholeinfirst;
	        int   em_isholeinsecond;
                float em_rookEnergy;
		float em_fracRook;
		float em_swissCross;
		float em_scRatio;
                int   em_hasBadSrpFlag; 
		float em_seedEnergy;
		float em_seedChi2;
		float em_seedTime;
                int   em_seedFlag;
                int   em_seedSrpFlag; 
		int   em_seedIphi;
		int   em_seedIeta;
		int   em_seedIx;
		int   em_seedIy; 
		int   em_seedZside;
		// energy after each energy corrections 
	        float emCorrEta_e;
                float emCorrEta_et;
	        float emCorrBR1_e;
	        float emCorrBR1_et;
	        float emCorrBR1Full_e;
	        float emCorrBR1Full_et;
		//photon informations
		int   em_isPhoton;
	        float pho_e;
        	float pho_et;
	        float pho_phi;
        	float pho_eta;
	        float pho_theta;
	        float pho_r9;
		float pho_HoE;
		// converted photons informations 
                int   pho_isConverted;
                int   pho_nTracks;
                float pho_EoverP;
                float pho_Rconv;
                float pho_Zconv;
                float pho_Xconv;
                float pho_Yconv;
                //electron informations
                int   em_isElectron;
                float ele_e;
                float ele_et;
                float ele_phi;
                float ele_eta;
                float ele_theta;
                float ele_charge;
                float ele_mass;
		float ele_EoverP;
		// the Monte Carlo truth 
		int   em_isMatchWithMC;
		int   mc_PDGType;
		float mc_e;
		float mc_et;
		float mc_eta;
		float mc_phi;
		float mc_theta;
		int   mc_nbMatch;
		// the MC truth for photon
		int   mc_isPhoton;
		int   mc_isConverted;
		float mc_convEt;
		float mc_convR;
		float mc_convX;
		float mc_convY;
		float mc_convZ;
		int   mc_Nconv;

      };
      tree_structure_ tree_;

      struct basicClusters_structure_ {
                //references to the Event
                int eventRef;
                int runNum;
                int bx;
                int orbite;
                int triggerType;
                int lumiBlock;
                // technical L1 trigger
                int techTrigger0;
                int techTrigger40;
                int techTrigger41;
                int techTrigger36;
                int techTrigger37;
                int techTrigger38;
                int techTrigger39;	
		// BC informations
                int   bc_isInCrack;
                int   bc_barrelOrEndcap;
                float bc_e;
                float bc_et;
                float bc_phi;
                float bc_eta;  
		// SC informations
		int   bc_isInSC;
		int   bc_refSC;
		float bc_fraction;
		// clusterShape info
		float bc_r9;
		int   bc_nbcrystal;
      };	
      basicClusters_structure_ treeBC_; 	

      struct event_structure{
                //references to the Event
                int eventRef;
                int runNum;
                int bx;
                int orbite;
                int triggerType;
                int lumiBlock;
              // technical L1 trigger
                int techTrigger0;
                int techTrigger40;
                int techTrigger41;
                int techTrigger36;
                int techTrigger37;
                int techTrigger38;
                int techTrigger39;
	     // BC SC containement
		int nbSuperClusterBarrel;
		int nbSuperClusterEndcapP;
		int nbSuperClusterEndcapM;
		int nbBasicClusterBarrel;
		int nbBasicClusterEndcapP;
		int nbBasicClusterEndcapM;
             // Selective ReadOut State
                int nbSrpTTEB;
		int nbZsTTEB;
		int nbFlag0TTEB; 
                int nbFlag2TTEB;
                int nbFlag4TTEB;
                int nbFlag5TTEB;
                int nbFlag6TTEB;
                int nbFlag7TTEB;
                int nbSrpTTEEM;
                int nbZsTTEEM;
		int nbFlag0TTEEM;	
                int nbFlag2TTEEM;
                int nbFlag4TTEEM;
                int nbFlag5TTEEM;
                int nbFlag6TTEEM;
                int nbFlag7TTEEM;
                int nbSrpTTEEP;
                int nbZsTTEEP;
		int nbFlag0TTEEP;
                int nbFlag2TTEEP;
                int nbFlag4TTEEP;
                int nbFlag5TTEEP;
                int nbFlag6TTEEP;
                int nbFlag7TTEEP;
		//nb RH EB
		int nbRH_EB;
		int nbRH_EE;


      };		
      event_structure treeEvent_;

      struct recHits_structure{
                //references to the Event
                int eventRef;
                int runNum;
                int bx;
                int orbite;
                int triggerType;
                int lumiBlock;
              // technical L1 trigger
                int techTrigger0;
                int techTrigger40;
                int techTrigger41;
                int techTrigger36;
                int techTrigger37;
                int techTrigger38;
                int techTrigger39;
		int techTrigger42;
		int techTrigger43;
              // RecHits Infos
		int rh_barrelOrEndcap;
		float rh_phi;
		float rh_eta; 
                int rh_iPhi;
                int rh_iEta;
                int rh_iSm;
                int rh_iSmPhi;
                int rh_iSmEta;
                int rh_iX;
                int rh_iY;
                int rh_zSide;
                float rh_energy;
		float rh_uncalibEnergy;
                float rh_chi2;
		float rh_time;
                int rh_flag;
                int rh_srpFlag;
                int rh_isCluster;
                int rh_isSuperCluster;
		float rh_eHfNeg;
		float rh_eHfPos;
		float rh_eHfNegTime;
		float rh_eHfPosTime;
		int   rh_eHfNcounts;
		int   rh_eHfPcounts; 
      }; 
      recHits_structure treeRecHits_;	

      struct trigger_structure{
                //references to the Event
                int eventRef;
                int runNum;
                int bx;
                int orbite;
                int triggerType;
                int lumiBlock;
              // technical L1 trigger
                int techTrigger0;
                int techTrigger40;
                int techTrigger41;
                int techTrigger36;
                int techTrigger37;
                int techTrigger38;
                int techTrigger39;
		// TT infos
 		int TT_barrelOrEndcap;
		int TT_iphi;
		int TT_ieta;
		int TT_ix;
		int TT_iy;
		int TT_zside;
		int TT_flag;	
      };
      trigger_structure treeTrigger_;

      struct HF_structure{
                //references to the Event
                int eventRef;
                int runNum;
                int bx;
                int orbite;
                int triggerType;
                int lumiBlock;
              // technical L1 trigger
                int techTrigger0;
                int techTrigger40;
                int techTrigger41;
                int techTrigger36;
                int techTrigger37;
                int techTrigger38;
                int techTrigger39;
	      // info on the HF recHit 	
		float hf_eta;
		float hf_phi;
		int   hf_ieta;
		int   hf_iphi;
                float hf_energy;
                float hf_time;
		int   hf_depth;
		float hf_alphaRatio;
		float hf_alphaRatioTimed;	
	};
	HF_structure treeHF_;
		 
};
