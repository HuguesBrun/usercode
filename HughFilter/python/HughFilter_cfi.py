import FWCore.ParameterSet.Config as cms
HughFilter = cms.EDFilter("HughFilter",
   endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
   barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
   superClusterCollection = cms.string(''),
   superClusterProducer = cms.string('correctedHybridSuperClusters'),
   superClusterCollectionEndcap = cms.string(''),
   superClusterProducerEndcap = cms.string('correctedMulti5x5SuperClustersWithPreshower'),
   L1triggerResults = cms.InputTag('gtDigis'),
   ScetThreshold = cms.double(0),
   NbXtalThreshold = cms.int32(0), 
   nEvent = cms.int32(-1),
   runNumber = cms.int32(124120),
   eventNumbers = cms.vint32(123321), 	
   keepOnlyTwenty = cms.bool(False)
)
