import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
#process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FakeConditions_cff")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/hbrun/theRECOfile.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.singleEModule = cms.EDAnalyzer("RecoPhotonEnergyScaleAnalyzer",
    outputFile = cms.string('Two_gamma2.root'),
    endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    correctedSuperClusterCollection = cms.string(''),
    superClusterProducer = cms.string('correctedHybridSuperClusters'),
    hitCollection = cms.string('EcalRecHitsEB'),
    correctedSuperClusterProducer = cms.string('correctedHybridSuperClusters'),
    correctedPhotonCollection = cms.string('photons'),
    hepMCLabel = cms.string('generator'),
    superClusterCollection = cms.string('')
)

process.p = cms.Path(process.singleEModule)


