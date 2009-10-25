import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
#process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FakeConditions_cff")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/hbrun/theRECOfile.root')
)

process.singleEModule = cms.EDAnalyzer("RecoPhotonEndcapESAnalyzer",
    outputFile = cms.string('DiPhotonEE.root'),
    endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    correctedSuperClusterCollection = cms.string(''),
    superClusterProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower'),
    correctedSuperClusterProducer = cms.string('correctedMulti5x5SuperClustersWithPreshower'),
    correctedPhotonCollection = cms.string('photons'),
    hepMCLabel = cms.string('generator'),
    superClusterCollection = cms.string('')
)

process.p = cms.Path(process.singleEModule)


