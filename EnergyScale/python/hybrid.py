import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")
process.load("Configuration.StandardSequences.Geometry_cff")


process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/30X/RECO/test_new.root')
)

process.singleEModule = cms.EDAnalyzer("EnergyScaleAnalyzer",
    recalculateWidths = cms.bool(False),
    outputFile = cms.string('hybrid_300.root'),
    endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),                                       
    hitProducer = cms.string('ecalRecHit'),
    correctedSuperClusterCollection = cms.string(''),
    superClusterProducer = cms.string('correctedHybridSuperClusters'),
    hitCollection = cms.string('EcalRecHitsEB'),
    correctedSuperClusterProducer = cms.string('correctedHybridSuperClusters'),
    hepMCLabel = cms.string('source'),
    superClusterCollection = cms.string('')
)

process.p = cms.Path(process.singleEModule)


