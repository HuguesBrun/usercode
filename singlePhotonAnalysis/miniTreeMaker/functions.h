#include "TMath.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7_new/src/Morgan/IpnTreeProducer/interface/TRootRun.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7_new/src/Morgan/IpnTreeProducer/interface/TRootPhoton.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7_new/src/Morgan/IpnTreeProducer/interface/TRootSuperCluster.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7_new/src/Morgan/IpnTreeProducer/interface/TRootElectron.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7_new/src/Morgan/IpnTreeProducer/interface/TRootJet.h"
#include "/sps/cms/hbrun/CMSSW_3_9_7_new/src/Morgan/IpnTreeProducer/interface/TRootMCParticle.h"


double DeltaR(double , double , double , double );

//int* InitializeHLTinfo(TRootRun* runInfos, string* ListWantedHLTnames, int nPathWanted);

double GetClosestJetEMFraction(TRootPhoton* myphoton, TClonesArray* jets, double JetPt);

double GetDeltaRClosestJet(TRootPhoton* myphoton, TClonesArray* jets, double JetPt);

double GetPtOverJetPt(TRootPhoton* myphoton, TClonesArray* jets, double JetPt);

double GetTransverseMomentumToJetDirection(TRootPhoton *myPhoton, TClonesArray* jets, double JetPt);

void doGenInfo(TRootPhoton* myphoton, TClonesArray* mcParticles, Int_t* pho_GenId, Int_t* pho_MotherId, Int_t* pho_isGenElectron, Int_t* pho_isPromptGenPho, Int_t* pho_isFromQuarkGen, Int_t* pho_isPi0Gen, Int_t* pho_isEtaGen, Int_t* pho_isRhoGen, Int_t* pho_isOmegaGen, Float_t* pho_PromptGenIsoEnergyStatus1, Float_t* pho_PromptGenIsoEnergyStatus2, double dRcone);

void doGenInfo(TRootPhoton* myphoton, TClonesArray* mcParticles, Int_t* pho_GenId, Int_t* pho_MotherId, Int_t* pho_isGenElectron, Int_t* pho_isPromptGenPho, Int_t* pho_isFromQuarkGen, Int_t* pho_isPi0Gen, Int_t* pho_isEtaGen, Int_t* pho_isRhoGen, Int_t* pho_isOmegaGen, Float_t* pho_PromptGenIsoEnergyStatus1, Float_t* pho_PromptGenIsoEnergyStatus2, float* pho_trueE, float* pho_truePx, float* pho_truePy, float* pho_truePz, float* pho_trueEta, float* pho_truePhi,double dRcone);
void matchWithAnElectron(TRootPhoton *myPhoton, TClonesArray *electrons, int *isAlsoaRecoElectron, float *pho_fBrem, float *pho_momentumCorrected, float *pho_d0, float *pho_tightEleId, float *pho_eleTrkIso, float *pho_eleEcalIso, float *pho_eleHcalIso, float *pho_eleDeltaPhiIn, float *pho_eleDeltaEtaIn, float *pho_eleHoE, float *pho_eleSigmaIetaIeta, int *pho_eleMissHits, float *pho_eleDistConvPartner, float *pho_eleDcotConvPartner);
