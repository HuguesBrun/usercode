#include "functions.h"
#include "functions.cc"


#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"

#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootBardak.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootBeamSpot.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootBeamStatus.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootCluster.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootDummyEvent.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootEcalRecHit.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootElectron.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootEvent.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootJet.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootMCParticle.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootMCPhoton.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootMET.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootMuon.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootParticle.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootPhoton.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootRun.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootSignalEvent.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootSuperCluster.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootTopTop.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootTrack.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootVertex.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7/src/Morgan/IpnTreeProducer/interface/TRootHLTObject.h"


TFile *myFile;// = new TFile("theMiniTree.root","RECREATE");
TTree *myTree_;
TChain *inputEventTree = new TChain("eventTree");
TChain *inputRunTree = new TChain("runTree");

//string ListWantedHLTnames[7] = {"HLT_DoublePhoton5_CEP_L1R","HLT_DoublePhoton17_L1R","HLT_Photon20_Cleaned_L1R","HLT_Photon30_Cleaned_L1R","HLT_Photon50_NoHE_Cleaned_L1R","HLT_DoublePhoton10_L1R","HLT_DoublePhoton15_L1R"};
string ListWantedHLTnames[12] = {"HLT_DoublePhoton22_L1R_v1","HLT_DoublePhoton17_SingleIsol_L1R_v1","HLT_Photon20_Cleaned_L1R","HLT_Photon30_Cleaned_L1R","HLT_Photon40_Isol_Cleaned_L1R_v1","HLT_DoublePhoton5_CEP_L1R_v3","HLT_Photon110_NoHE_Cleaned_L1R_v1","HLT_Photon17_Isol_SC17HE_L1R_v1","HLT_Photon22_SC22HE_L1R_v1","HLT_Photon40_CaloId_Cleaned_L1R_v1","HLT_Photon50_Cleaned_L1R_v1","HLT_Photon70_Cleaned_L1R_v1"};
int nbHlt = 12;

  bool doHLT;
  bool doHLTobject;
  bool doMC;
  bool doJetMC;
  bool doMETMC;
  bool doPDFInfo;
  bool doSignalMuMuGamma;
  bool doSignalTopTop;
  bool doPhotonConversionMC;
  bool doBeamSpot;
  bool doPrimaryVertex;
  bool doZeePrimaryVertex;
  bool doTrack; 
  bool doJet;
  bool doMuon;
  bool doElectron;
  bool doPhoton;
  bool doCluster;
  bool doPhotonConversion;
  bool doMET;
  bool doBardak;
  bool doPhotonVertexCorrection;
  bool doPhotonIsolation;    


                int pho_iEvent;
		int pho_EventNumber;
		int pho_RunNumber;
		int pho_lumiSection;

                int pho_isEB;
                int pho_isEE;
                int pho_isEEP;
                int pho_isEEM;
                int pho_isAlsoElectron;
                int pho_ElectronClassification;
                int pho_hasPixelSeed;
				int pho_eventProcessId;
				float pho_et;
                float pho_energy;
                float pho_px;
                float pho_py;
                float pho_pz;
                float pho_eta;
                float pho_phi;
                float pho_s4;
                float pho_s9;
                float pho_eMax;
                float pho_e2nd;
                float pho_r19;
                float pho_r9;
                float pho_S9overS9minusS1S2;
                float pho_cEE;
                float pho_cEP;
                float pho_cPP;
                float pho_phiwidth;
                float pho_etawidth;
                float pho_sigmaphi;
                float pho_sigmaeta;
                float pho_sigmaIetaIeta;
				float pho_sigmaEtaEta;
                float pho_hoe;
                float pho_IsoEcalRechit;
                float pho_IsoHcalRechit;
                float pho_IsoSolidTrkCone;
                float pho_IsoHollowTrkCone;
                float pho_IsoEcalRechit03;
                float pho_IsoHcalRechit03;
                float pho_IsoSolidTrkCone03;
                float pho_IsoHollowTrkCone03;
                float pho_esRatio;
                int pho_convNTracks;
                float pho_ptoverjetpt;
                float pho_DrSCclosest;
                float pho_DrTrkclosest;
                float pho_ptoverjetpt_pt2;
                float pho_DrTrkclosest_pt2;
                float pho_ptoverjetpt_pt5;
                float pho_DrTrkclosest_pt5;
                float pho_ptoverjetpt_pt10;
                float pho_DrTrkclosest_pt10;
                float pho_ptoverjetpt_pt15;
                float pho_ptoverjetpt_PFlow;
                float pho_ptoverjetpt_PFlow_pt2;
                float pho_ptoverjetpt_PFlow_pt5;
                float pho_ptoverjetpt_PFlow_pt10;
                float pho_ptoverjetpt_sisCone;
                float pho_ptoverjetpt_sisCone_pt2;
                float pho_ptoverjetpt_sisCone_pt5;
                float pho_ptoverjetpt_sisCone_pt10;
				float pho_transverseMomentumToJetDirection;
				float pho_transverseMomentumToJetDirection_pt2;
				float pho_transverseMomentumToJetDirection_pt5;
				float pho_transverseMomentumToJetDirection_pt10;
				float pho_transverseMomentumToJetDirection_pt20;
				float pho_transverseToJetRatio;
				float pho_transverseToJetRatio_pt2;
				float pho_transverseToJetRatio_pt5;
				float pho_transverseToJetRatio_pt10;
				float pho_transverseToJetRatio_pt20;
				// other collection matching
				int   pho_isAlsoRecoAsElectron;
				float pho_fBrem;
				float pho_momentumCorrected;
				float pho_d0;
				int   pho_isAlsoRecoAsJet;		
		float pho_matchWithEle;
		float pho_tightEleId;
		float pho_eleTrkIso;
		float pho_eleEcalIso;
		float pho_eleHcalIso;
		float pho_eleDeltaPhiIn;
		float pho_eleDeltaEtaIn;
		float pho_eleHoE;
		float pho_eleSigmaIetaIeta;
		int   pho_eleMissHits;
		float pho_eleDistConvPartner;
		float pho_eleDcotConvPartner;
                float pho_jetEMfraction;
                float pho_DrJetClosest;
                int pho_isLoose;
                int pho_isTight;
                int pho_isEBGap;
                int pho_isEEGap;
                int pho_GenId;
                int pho_MotherId;
                int pho_isPromptGenPho;
                int pho_isFromQuarkGen;
                int pho_isPi0Gen;
                int pho_isEtaGen;
                int pho_isRhoGen;
                int pho_isOmegaGen;
                int pho_isGenElectron;
                int pho_eventPassHLT_Photon10_L1R;
                int pho_eventPassHLT_Photon15_L1R;
                int pho_eventPassHLT_DoublePhoton10_L1R;
                float pho_eventPtHat;
                int pho_nVertex;
                float pho_PromptGenIsoEnergyStatus1_cone02;
                float pho_PromptGenIsoEnergyStatus2_cone02;
                float pho_PromptGenIsoEnergyStatus1_cone03;
                float pho_PromptGenIsoEnergyStatus2_cone03;
                float pho_PromptGenIsoEnergyStatus1_cone035;
                float pho_PromptGenIsoEnergyStatus2_cone035;
                float pho_PromptGenIsoEnergyStatus1_cone04;
                float pho_PromptGenIsoEnergyStatus2_cone04;
                int pho_seedSeverity;
                int pho_recoFlag;
		float pho_seedTime;
		float pho_seedEnergy;
		int pho_HLT_bit0;
		int pho_HLT_bit1;
		int pho_HLT_bit2;
		int pho_HLT_bit3;
		int pho_HLT_bit4;
		int pho_HLT_bit5;
		int pho_HLT_bit6;
		int pho_HLT_bit7;
		int pho_HLT_bit8;
		int pho_HLT_bit9;
		int pho_HLT_bit10;
		int pho_HLT_bit11;
		// sc infos
		float pho_SCeta;
		float pho_SCphi;
		float pho_SCEtraw;
		float pho_SCEt;
		float pho_SCr9;
		float pho_SCbr;
		int   pho_SCnbBC;
		int   pho_SCnXtal;
		int  isAspike; // if we want do ID spikes ;)	
		float pho_etaLAT;
		float pho_phiLAT;
		float pho_LAT;
		float pho_Zernike20;
		float pho_Zernike42;
		float pho_secondMomentMaj;
		float pho_secondMomentMin;
		float pho_secondMomentAlpha;
		float pho_ESratio;	
		float pho_trueE;
		float pho_truePx;
		float pho_truePy;
		float pho_truePz;
		float pho_trueEta;
		float pho_truePhi;
		// HLT object
        int pho_isMatchingWithHLTObject;

TBranch* event_br = 0;
TBranch* run_br = 0;
TBranch* mcParticles_br = 0;
TBranch* genJets_br = 0;
TBranch* genMETs_br = 0;
TBranch* mcSignalMuMuGamma_br = 0;
TBranch* mcTopTopEvent_br = 0;
TBranch* mcPhotons_br = 0;
TBranch* beamSpot_br = 0;
TBranch* vertices_br = 0;
TBranch* zeeVertices_br = 0;
TBranch* tracks_br = 0;
TBranch* jets_br = 0;
TBranch* pflowjets_br = 0;
TBranch* sisconejets_br = 0;
TBranch* muons_br = 0;
TBranch* electrons_br = 0;
TBranch* photons_br = 0;
TBranch* clusters_br = 0;
TBranch* superClusters_br = 0;
TBranch* conversions_br = 0;
TBranch* met_br = 0;
TBranch* bardak_br = 0;
TBranch* HLTObjects_br = 0;

TRootEvent* event = 0;
TRootRun* runInfos = 0;
TClonesArray* mcParticles; // = new TClonesArray("TRootMCParticle", 0);
TClonesArray* genJets; // = new TClonesArray("TRootParticle", 0);
TClonesArray* genMETs; // = new TClonesArray("TRootParticle", 0);
TRootSignalEvent* mcMuMuGammaEvent = 0;
TRootSignalEvent* mcTopTopEvent = 0;
TClonesArray* mcPhotons; // = new TClonesArray("TRootMCPhoton", 0);
TRootBeamSpot* beamSpot = 0;
TClonesArray* vertices; // = new TClonesArray("TRootVertex", 0);
TClonesArray* zeeVertices; // = new TClonesArray("TRootVertex", 0);
TClonesArray* tracks; // = new TClonesArray("TRootTrack", 0);
TClonesArray* jets; // = new TClonesArray("TRootJet", 0);
TClonesArray* pflowjets; // = new TClonesArray("TRootJet", 0);
TClonesArray* sisconejets; // = new TClonesArray("TRootJet", 0);
TClonesArray* muons; /// = new TClonesArray("TRootMuon", 0);
TClonesArray* electrons; // = new TClonesArray("TRootElectron", 0);
TClonesArray* photons; // = new TClonesArray("TRootPhoton", 0);
TClonesArray* clusters; // = new TClonesArray("TRootCluster", 0);
TClonesArray* superClusters; // = new TClonesArray("TRootSuperCluster", 0);
TClonesArray* conversionTracks; // = new TClonesArray("TRootTrack", 0);
TClonesArray* met; // = new TClonesArray("TRootMET", 0);
TRootBardak* bardak = 0;
TRootSuperCluster* mysc;
TClonesArray* HLTObjects;
