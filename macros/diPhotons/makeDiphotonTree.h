#include "functions.h"
#include "functions.C"

#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"
#include <TLorentzVector.h>

TFile *myFile;// = new TFile("theMiniTree.root","RECREATE");
TTree *myTree_;
TChain *chain = new TChain("photons");

  bool doMC;
  int NbEntries;

// MiniTtree d'entr�e
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
                int pho_isAfterCut1;
                int pho_isAfterCut2;
                int pho_isAfterCut3;
                int pho_isAfterCut4;
                int pho_isAfterCut5;
                int pho_isAfterCut6;
                int pho_isAfterCut7;
                int pho_isAfterCut8;
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
                int pho_HLT30;
                int pho_HLT50;
                int pho_HLTdiPho15;
                int pho_HLTdiPho20;
                // sc infos
                float pho_SCeta;
                float pho_SCphi;
                float pho_SCEtraw;
                float pho_SCEt;
                float pho_SCr9;
                float pho_SCbr;
                int   pho_SCnbBC;
                int   pho_SCnXtal;

// miniTree de sortie

int event_number;    
int event_realNumber;
int event_runNumber;                                    
int event_LumiSection; 
int event_weightLO;                                       
int event_weightNLO;                                     
int event_eventPtHat; 
int event_nVertex;
int dipho_cat4;                                           
int dipho_cat4bis;                                        
int dipho_cat6;                                           
int dipho_cat12;                                          
float dipho_mgg;                                
float dipho_qt;                                           
float dipho_ql;                                           
float dipho_deltaR;                                       
float dipho_costhetastar;                                 
float dipho_eta;                                          
float dipho_etastar;                                      
float dipho_rapidity;                                     
float dipho_deltaphi;                                     
float dipho_costhetastar_parton;                          
float dipho_tanhrapiditydiff;                             
float pholead_pt;                                         
float photrail_pt;                                        
float pho1_pt;                                            
float pho2_pt;                                            
float pholead_eta;                                        
float photrail_eta;                                       
float pholead_r9;                                         
float photrail_r9;                                        
float pholead_cPP;                                        
float photrail_cPP;                                       
float pholead_sigmaphi;                                   
float photrail_sigmaphi;                                  
float pholead_S9overS9minusS1S2;                          
float photrail_S9overS9minusS1S2;                         
float pholead_etawidth;                                   
float photrail_etawidth;                                  
float pholead_sigieta;                                    
float photrail_sigieta;                                   
float pholead_ptoverjetpt;                                
float photrail_ptoverjetpt;                               
float pholead_DrScClosest;                                
float photrail_DrScClosest;                               
float pholead_DrTrkClosest;                               
float photrail_DrTrkClosest;                              
float pholead_MLPoutputClusterShape;                      
float photrail_MLPoutputClusterShape;                     
float pholead_MLPoutputEnv;                               
float photrail_MLPoutputEnv;                              
float pholead_MLPoutputClusterShapeEnv;                   
float photrail_MLPoutputClusterShapeEnv;                  
float phoNNmin_MLPoutputClusterShape;                     
float phoNNmax_MLPoutputClusterShape;                     
float phoNNmin_MLPoutputEnv;                              
float phoNNmax_MLPoutputEnv;                              
float phoNNmin_MLPoutputClusterShapeEnv;                  
float phoNNmax_MLPoutputClusterShapeEnv;                  
float pholead_MLPoutputNNisol;                            
float photrail_MLPoutputNNisol;                           
float pho1_MLPoutputNNisol;                               
float pho2_MLPoutputNNisol;                               
float phoNNmin_MLPoutputNNisol;                           
float phoNNmax_MLPoutputNNisol;
float pholead_HcalIso;
float photrail_HcalIso;
float pholead_EcalIso;
float photrail_EcalIso;
float pholead_TrackerIso;
float photrail_TrackerIso;
float pholead_HcalIsoPerso;
float photrail_HcalIsoPerso;
float pholead_EcalIsoPerso;
float photrail_EcalIsoPerso;
float pholead_TrackerIsoPerso;
float photrail_TrackerIsoPerso;
float pholead_hoe;
float photrail_hoe;
float pholead_DrJetClosest;
float photrail_DrJetClosest;
float pholead_jetEMfraction;
float photrail_jetEMfraction;
int pholead_isPromptGenPho;
int photrail_isPromptGenPho;
int pholead_isFromQuarkGen;
int photrail_isFromQuarkGen;
int pholead_HasPixSeed;
int photrail_HasPixSeed;
int pholead_seedSeverity;
int photrail_seedSeverity;
int pholead_recoFlag;
int photrail_recoFlag;
int pholead_isEB;
int photrail_isEB;
int pholead_isEE;
int photrail_isEE;
int dipho_eventPassHLT_Photon30_L1R;
int dipho_eventPassHLT_Photon50_L1R;
int dipho_eventPassHLT_DoublePhoton10_L1R;
int dipho_eventPassHLT_DoublePhoton15_L1R;
float dipho_Kfactor;


  TLorentzVector P_loc[100] ;
  double Eta_loc[100] ;
  double R9_loc[100] ;
  double cPP_loc[100] ;
  double EtaWidth_loc[100] ;
  double SigIeta_loc[100] ;
  double sigmaphi_loc[100] ;
  double S9overS9minusS1S2_loc[100] ;
  double ptoverjetpt_loc[100] ;
  double closestSC_dR_loc[100] ;
  double dR_SCtrkclosest_loc[100] ;

  double MLPoutputClusterShape_loc[100] ;
  double MLPoutputEnv_loc[100] ;
  double MLPoutputClusterShapeEnv_loc[100] ;
  double MLPoutputNNisol_loc[100] ;

  double TrackerIso_loc[100] ;
  double EcalIso_loc[100] ;
  double HcalIso_loc[100] ;
  double TrackerIsoPerso_loc[100] ;
  double EcalIsoPerso_loc[100] ;
  double HcalIsoPerso_loc[100] ;
  double HoE_loc[100] ;
  int HasPixSeed_loc[100] ;
  int SeedSeverity_loc[100] ;
  int RecoFlag_loc[100] ;

  double DrJetClosest_loc[100] ;
  double JetEMfraction_loc[100] ;

  int isPromptGenPho_loc[100] ;
  int isFromQuarkGen_loc[100] ;

  int passHLT_Photon10_loc[100] ;
  int passHLT_Photon15_loc[100] ;
  int passHLT_DoublePhoton10_loc[100] ;

  int pho_HLT30_loc[100] ;
  int pho_HLT50_loc[100] ;
  int pho_HLTdiPho15_loc[100] ;
  int pho_HLTdiPho20_loc[100] ;

  double Kfactor_loc[100] ;

  float pho_sigmaIetaIeta_loc[100];
  int pho_isEB_loc[100];
  int pho_isEE_loc[100];

