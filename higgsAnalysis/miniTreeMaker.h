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
#include <utility>

#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootBardak.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootBeamSpot.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootBeamStatus.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootCluster.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootDummyEvent.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootEcalRecHit.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootElectron.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootEvent.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootJet.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMCParticle.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMCPhoton.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMET.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMuon.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootParticle.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootPhoton.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootRun.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootSignalEvent.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootSuperCluster.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootTopTop.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootTrack.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootVertex.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootHLTObject.h"


TFile *myFile;// = new TFile("theMiniTree.root","RECREATE");
TTree *myTree_;
TChain *inputEventTree = new TChain("eventTree");
TChain *inputRunTree = new TChain("runTree");

//string ListWantedHLTnames[13] = {"HLT_DoublePhoton33_v1","HLT_Photon125_NoSpikeFilter_v1","HLT_Photon20_R9Id_Photon18_R9Id_v1","HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1","HLT_Photon26_CaloIdL_IsoVL_Photon18_v1","HLT_Photon26_IsoVL_Photon18_IsoVL_v1","HLT_Photon26_IsoVL_Photon18_v1","HLT_Photon26_Photon18_v1","HLT_Photon30_CaloIdVL_IsoL_v1","HLT_Photon30_CaloIdVL_v1","HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1","HLT_Photon75_CaloIdVL_IsoL_v1","HLT_Photon75_CaloIdVL_v1"};

//string ListWantedHLTnames[18] = {"HLT_DoublePhoton33_v2","HLT_Photon125_NoSpikeFilter_v2","HLT_Photon20_CaloIdVL_IsoL_v1","HLT_Photon20_R9Id_Photon18_R9Id_v2","HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2","HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v1","HLT_Photon26_CaloIdL_IsoVL_Photon18_v2","HLT_Photon26_IsoVL_Photon18_IsoVL_v2","HLT_Photon26_IsoVL_Photon18_v2","HLT_Photon26_Photon18_v2","HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v1","HLT_Photon30_CaloIdVL_IsoL_v2","HLT_Photon30_CaloIdVL_v2","HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2","HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1","HLT_Photon50_CaloIdVL_IsoL_v1","HLT_Photon75_CaloIdVL_IsoL_v2","HLT_Photon75_CaloIdVL_v2"};

//string ListWantedHLTnames[18] = {"HLT_DoublePhoton33_v3","HLT_Photon125_NoSpikeFilter_v3","HLT_Photon20_CaloIdVL_IsoL_v2","HLT_Photon20_R9Id_Photon18_R9Id_v3","HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3","HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v2","HLT_Photon26_CaloIdL_IsoVL_Photon18_v3","HLT_Photon26_IsoVL_Photon18_IsoVL_v3","HLT_Photon26_IsoVL_Photon18_v3","HLT_Photon26_Photon18_v3","HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v2","HLT_Photon30_CaloIdVL_IsoL_v3","HLT_Photon30_CaloIdVL_v3","HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3","HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2","HLT_Photon50_CaloIdVL_IsoL_v2","HLT_Photon75_CaloIdVL_IsoL_v3","HLT_Photon75_CaloIdVL_v3"};

/*string ListWantedHLTnames[34] = { // Run2011 1e33 v1.3 
"HLT_DoubleEle33_CaloIdL_v1",
"HLT_DoubleEle33_v1",
"HLT_DoublePhoton33_HEVT_v1",
"HLT_DoublePhoton33_v4",
"HLT_DoublePhoton40_MR150_v1",
"HLT_DoublePhoton40_R014_MR150_v1",
"HLT_DoublePhoton50_v1",
"HLT_DoublePhoton5_IsoVL_CEP_v3",
"HLT_DoublePhoton60_v1",
"HLT_Photon125_v1",
"HLT_Photon200_NoHE_v1",
"HLT_Photon20_CaloIdVL_IsoL_v3",
"HLT_Photon20_R9Id_Photon18_R9Id_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v4",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v4",
"HLT_Photon26_IsoVL_Photon18_v4",
"HLT_Photon26_Photon18_v4",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3",
"HLT_Photon26_R9Id_Photon18_R9Id_v1",
"HLT_Photon30_CaloIdVL_IsoL_v4",
"HLT_Photon30_CaloIdVL_v4",
"HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v1",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3",
"HLT_Photon36_IsoVL_Photon22_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1",
"HLT_Photon50_CaloIdVL_IsoL_v3",
"HLT_Photon50_CaloIdVL_v1",
"HLT_Photon75_CaloIdVL_IsoL_v4",
"HLT_Photon75_CaloIdVL_v4",
"HLT_Photon90_CaloIdVL_IsoL_v1",
"HLT_Photon90_CaloIdVL_IsoL_v1"
};*/


/*string ListWantedHLTnames[38] = { // Run2011 1e33 v2.3
"HLT_DoubleEle33_CaloIdL_v2",
"HLT_DoubleEle33_v2",
"HLT_DoublePhoton33_HEVT_v2",
"HLT_DoublePhoton33_v5",
"HLT_DoublePhoton40_MR150_v3",
"HLT_DoublePhoton40_R014_MR150_v3",
"HLT_DoublePhoton50_v2",
"HLT_DoublePhoton5_IsoVL_CEP_v4",
"HLT_DoublePhoton60_v2",
"HLT_Photon125_v2",
"HLT_Photon200_NoHE_v2",
"HLT_Photon20_CaloIdVL_IsoL_v4",
"HLT_Photon20_R9Id_Photon18_R9Id_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v5",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v5",
"HLT_Photon26_IsoVL_Photon18_v5",
"HLT_Photon26_Photon18_v5",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4",
"HLT_Photon26_R9Id_Photon18_R9Id_v2",
"HLT_Photon30_CaloIdVL_IsoL_v5",
"HLT_Photon30_CaloIdVL_v5",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v2",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4",
"HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1",
"HLT_Photon36_IsoVL_Photon22_v2",
"HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_R9Id_Photon22_R9Id_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2",
"HLT_Photon50_CaloIdVL_IsoL_v4",
"HLT_Photon50_CaloIdVL_v2",
"HLT_Photon75_CaloIdVL_IsoL_v5",
"HLT_Photon75_CaloIdVL_v5",
"HLT_Photon90_CaloIdVL_IsoL_v2",
"HLT_Photon90_CaloIdVL_v2"
};*/


/*string ListWantedHLTnames[38] = { // Run2011 1e33 v2.4
"HLT_DoubleEle33_CaloIdL_v2",
"HLT_DoubleEle33_v2",
"HLT_DoublePhoton33_HEVT_v2",
"HLT_DoublePhoton33_v5",
"HLT_DoublePhoton40_MR150_v3",
"HLT_DoublePhoton40_R014_MR150_v3",
"HLT_DoublePhoton50_v2",
"HLT_DoublePhoton5_IsoVL_CEP_v4",
"HLT_DoublePhoton60_v2",
"HLT_Photon125_v2",
"HLT_Photon200_NoHE_v2",
"HLT_Photon20_CaloIdVL_IsoL_v4",
"HLT_Photon20_R9Id_Photon18_R9Id_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v5",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v5",
"HLT_Photon26_IsoVL_Photon18_v5",
"HLT_Photon26_Photon18_v5",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4",
"HLT_Photon26_R9Id_Photon18_R9Id_v2",
"HLT_Photon30_CaloIdVL_IsoL_v5",
"HLT_Photon30_CaloIdVL_v5",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v2",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4",
"HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1",
"HLT_Photon36_IsoVL_Photon22_v2",
"HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1",
"HLT_Photon36_R9Id_Photon22_R9Id_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2",
"HLT_Photon50_CaloIdVL_IsoL_v4",
"HLT_Photon50_CaloIdVL_v2",
"HLT_Photon75_CaloIdVL_IsoL_v5",
"HLT_Photon75_CaloIdVL_v5",
"HLT_Photon90_CaloIdVL_IsoL_v2",
"HLT_Photon90_CaloIdVL_v2"
};*/

//string ListWantedHLTnames[12] = {"HLT_DoublePhoton22_L1R_v1","HLT_DoublePhoton17_SingleIsol_L1R_v1","HLT_Photon20_Cleaned_L1R","HLT_Photon30_Cleaned_L1R","HLT_Photon40_Isol_Cleaned_L1R_v1","HLT_DoublePhoton5_CEP_L1R_v3","HLT_Photon110_NoHE_Cleaned_L1R_v1","HLT_Photon17_Isol_SC17HE_L1R_v1","HLT_Photon22_SC22HE_L1R_v1","HLT_Photon40_CaloId_Cleaned_L1R_v1","HLT_Photon50_Cleaned_L1R_v1","HLT_Photon70_Cleaned_L1R_v1"};

string ListWantedHLTnames[44] = { // Monte Carlo 1e33
"HLT_DoubleEle33_CaloIdL_v3",
"HLT_DoubleEle33_v3",
"HLT_DoubleEle45_CaloIdL_v2",
"HLT_DoublePhoton33_HEVT_v3",
"HLT_DoublePhoton38_HEVT_v2",
"HLT_DoublePhoton40_MR150_v4",
"HLT_DoublePhoton40_R014_MR150_v4",	
"HLT_DoublePhoton5_IsoVL_CEP_v5",
"HLT_DoublePhoton60_v3",
"HLT_DoublePhoton80_v1",
"HLT_Photon135_v1",
"HLT_Photon200_NoHE_v3",
"HLT_Photon20_CaloIdVL_IsoL_v5",
"HLT_Photon20_R9Id_Photon18_R9Id_v6",
"HLT_Photon225_NoHE_v1",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v6",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v5",
"HLT_Photon26_CaloIdL_IsoVL_Photon18_v6",
"HLT_Photon26_IsoVL_Photon18_IsoVL_v6",
"HLT_Photon26_IsoVL_Photon18_v6",
"HLT_Photon26_Photon18_v6",
"HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v5",	
"HLT_Photon26_R9Id_Photon18_R9Id_v5",
"HLT_Photon30_CaloIdVL_IsoL_v6",
"HLT_Photon30_CaloIdVL_v6",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v2",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v2",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id_v1",
"HLT_Photon36_CaloIdL_IsoVL_Photon22_v3",
"HLT_Photon36_CaloIdL_Photon22_CaloIdL_v5",
"HLT_Photon36_CaloIdVL_Photon22_CaloIdVL_v1",
"HLT_Photon36_IsoVL_Photon22_v3",
"HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v2",
"HLT_Photon36_R9Id_Photon22_R9Id_v2",
"HLT_Photon400_v1",
"HLT_Photon40_CaloIdL_Photon28_CaloIdL_v3",
"HLT_Photon44_CaloIdL_Photon34_CaloIdL_v1",
"HLT_Photon48_CaloIdL_Photon38_CaloIdL_v1",
"HLT_Photon50_CaloIdVL_IsoL_v5",
"HLT_Photon50_CaloIdVL_v3",
"HLT_Photon75_CaloIdVL_IsoL_v6",	
"HLT_Photon75_CaloIdVL_v6",
"HLT_Photon90_CaloIdVL_IsoL_v3",
"HLT_Photon90_CaloIdVL_v3"
};
int nbHlt = 44; 







float secondPhotonCut = 20.0;
float IsoTrk = 6.0;
float IsoEcal = 4.0;
float IsoHcal = 4.0;


TString theHTLobject = "hltPhoton26IsoVLTrackIsoFilter";

  bool doHLT;
  bool doHLTobject;
  bool doMC;
  bool doJetMC;
  bool doMETMC;
  bool doPDFInfo;
  bool doDiphotons;
  bool doLeadingPhoton;
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
  bool doElectronConversionMC;
  bool doWorstIsolation; 

int dipho_HLT_bit0;
int dipho_HLT_bit1;
int dipho_HLT_bit2;
int dipho_HLT_bit3;
int dipho_HLT_bit4;
int dipho_HLT_bit5;
int dipho_HLT_bit6;
int dipho_HLT_bit7;
int dipho_HLT_bit8;
int dipho_HLT_bit9;
int dipho_HLT_bit10;
int dipho_HLT_bit11;
int dipho_HLT_bit12;
int dipho_HLT_bit13;
int dipho_HLT_bit14;
int dipho_HLT_bit15;
int dipho_HLT_bit16;
int dipho_HLT_bit17;
int dipho_HLT_bit18;
int dipho_HLT_bit19;
int dipho_HLT_bit20;
int dipho_HLT_bit21;
int dipho_HLT_bit22;
int dipho_HLT_bit23;
int dipho_HLT_bit24;
int dipho_HLT_bit25;
int dipho_HLT_bit26;
int dipho_HLT_bit27;
int dipho_HLT_bit28;
int dipho_HLT_bit29;
int dipho_HLT_bit30;
int dipho_HLT_bit31;
int dipho_HLT_bit32;
int dipho_HLT_bit33;
int dipho_HLT_bit34;
int dipho_HLT_bit35;
int dipho_HLT_bit36;
int dipho_HLT_bit37;
int dipho_HLT_bit38;
int dipho_HLT_bit39;
int dipho_HLT_bit40;
int dipho_HLT_bit41;
int dipho_HLT_bit42;
int dipho_HLT_bit43;
int dipho_HLT_bit44;
int dipho_HLT_bit45;
// event infos

// miniTree de sortie

//here the event infos
int event_number;    
int event_runNumber;                                    
int event_LumiSection; 
float event_eventPtHat; 
int event_processId;
int event_nRecoVertex;
int event_nGenInTimeVertex;
int event_nGenOutOfTimeVertex;
int event_nPhotons;
float  event_rho;

// diphoton kine
float dipho_mgg;                                
float dipho_qt;                                           
float dipho_ql;                                           
float dipho_deltaR;                                       
float dipho_costhetastar;                                 
float dipho_eta;
float dipho_etastar;

// MC diphoton kine
float diphoMC_mgg;                                
float diphoMC_qt;                                           
float diphoMC_ql;                                           
float diphoMC_deltaR;                                       
float diphoMC_costhetastar;                                 
float diphoMC_eta;
float diphoMC_etastar;





//phoLead
//relatedMC
int pholead_isMatchingWithMC;
// kinevars
float pholead_pt;
float pholead_eta;
float pholead_etaSC;
// cluster shape vars
float pholead_r9;
float pholead_cPP;
float pholead_cEP;
float pholead_cEE;
float pholead_r19;
float pholead_SCEraw;
float pholead_eMax;
float pholead_e2x2;
float pholead_e5x5;
float pholead_ratioSeed;
float pholead_ratioS4;
float pholead_lambdaRatio;
float pholead_lamdbaDivCov;
float pholead_lambdaRatioRand;
float pholead_lamdbaDivCovRand;
float pholead_sigmaphi;
float pholead_secondMomentMaj;
float pholead_secondMomentMin;
float pholead_secondMomentAlpha;
float pholead_covAngle;
float pholead_covAngle2;
float pholead_S9overS9minusS1S2;
float pholead_etawidth;
float pholead_sigieta;
float pholead_SCbr;

//isolation variables
float pholead_HcalIso;
float pholead_EcalIso;
float pholead_TrackerIso;
float pholead_HcalIsodR03;
float pholead_EcalIsodR03;
float pholead_TrackerIsodR03;
float pholead_HcalIsoPerso;
float pholead_EcalIsoPerso;
float pholead_hoe;

//other variables
int pholead_HasPixSeed;
int pholead_isPromptGenPho;
int pholead_isFromQuarkGen;
int pholead_seedSeverity;
int pholead_recoFlag;
int pholead_isEB;
int pholead_isEE;
float pholead_NNshapeOutput;




//phoTrail
//relatedMC
int photrail_isMatchingWithMC;


TBranch* event_br = 0;
TBranch* run_br = 0;
TBranch* mcParticles_br = 0;
TBranch* genJets_br = 0;
TBranch* genMETs_br = 0;
TBranch* mcSignalMuMuGamma_br = 0;
TBranch* mcTopTopEvent_br = 0;
TBranch* mcPhotons_br = 0;
TBranch* mcElectrons_br = 0;
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
TClonesArray* mcElectrons;
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
