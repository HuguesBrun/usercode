import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FakeConditions_cff")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/21X/RECO/test.root')
)

process.singleEModule = cms.EDAnalyzer("EnergyScaleAnalyzer",
    recalculateWidths = cms.bool(True),
    outputFile = cms.string('test.root'),
    endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    correctedSuperClusterCollection = cms.string(''),
    superClusterProducer = cms.string('correctedHybridSuperClusters'),
    hitCollection = cms.string('EcalRecHitsEB'),
    correctedSuperClusterProducer = cms.string('correctedHybridSuperClusters'),
    hepMCLabel = cms.string('source'),
    superClusterCollection = cms.string('')
)

process.p = cms.Path(process.singleEModule)


