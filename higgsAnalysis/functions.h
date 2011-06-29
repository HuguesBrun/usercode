#ifndef THECICCUT_H
#define THECICCUT_H
#include "theCiCCut.h"
#endif
#include "TMath.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootRun.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootPhoton.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootSuperCluster.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootElectron.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootJet.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMCParticle.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMCPhoton.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootMCElectron.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootHLTObject.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootVertex.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootTrack.h"
#include "/sps/cms/hbrun/CMSSW_4_2_3/src/Morgan/IpnTreeProducer/interface/TRootBeamSpot.h"


double DeltaR(double , double , double , double );

double min(double, double);
double localtrackIsolation(TRootVertex, TClonesArray *, TRootPhoton, TRootBeamSpot);
double  calcWorstTrackIsolation(TClonesArray*, TClonesArray*, TRootPhoton, TRootBeamSpot);
double dRtoTrack(TRootPhoton , TClonesArray* );
bool photonPassingPreselection(TRootPhoton);
bool photonIsPassingCIC(TRootPhoton, TClonesArray*, TClonesArray*, TRootBeamSpot, TClonesArray*);
double CosThetaStar(TLorentzVector, TLorentzVector);
int  findGenParticle(TRootPhoton *, TClonesArray *, TRootParticle *);
