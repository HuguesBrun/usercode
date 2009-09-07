import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/21X/RECO/test_new.root')
)

process.singleEModule = cms.EDAnalyzer("EndcapES",
    outputFile = cms.string('multi_214.root'),
    hitProducer = cms.string('ecalRecHit'),
    correctedSuperClusterCollection = cms.string(''),
    superClusterProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower'),
    hitCollection = cms.string('EcalRecHitsEB'),
    correctedSuperClusterProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower'),
    hepMCLabel = cms.string('source'),
    superClusterCollection = cms.string('')
)

process.p = cms.Path(process.singleEModule)


