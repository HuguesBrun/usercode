import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("hugues.HughFilter.HughFilter_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR10_P_V4::All'
#process.GlobalTag.globaltag = 'GR_R_35X_V6::All'

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#file:/tmp/hbrun/theRECOfile.root'
        '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root'
)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.singleEModule = cms.EDAnalyzer('FirstDataAnalyzer',
   readRecHits					=	cms.bool(False),  
   readTTinfos					=	cms.bool(False), 
   readHFrecHits				=	cms.bool(False),
   readMCTruth					=	cms.bool(True),
   readPhotonMCTruth				=	cms.bool(True),
   deltaRMax					= 	cms.double(0.3),
   L1triggerResults 				= 	cms.InputTag('gtDigis'),
   HFrecHitsCollection				= 	cms.InputTag('hfreco'),
   barrelSrpFlagsCollection			= 	cms.string('ecalDigis'),
   endcapSrpFlagsCollection                     =       cms.string('ecalDigis'),
   barrelEcalHits   				= 	cms.InputTag('ecalRecHit','EcalRecHitsEB'),
   endcapEcalHits   				= 	cms.InputTag('ecalRecHit','EcalRecHitsEE'),
   barrelClusterProducer	      		=	cms.string('hybridSuperClusters'),
   barrelClusterCollection			=	cms.string('hybridBarrelBasicClusters'),
   endcapClusterProducer			=	cms.string('multi5x5BasicClusters'),
   endcapClusterCollection			=	cms.string('multi5x5EndcapBasicClusters'),
   barrelCorrectedSuperClusterCollection	= 	cms.string(''),
   barrelCorrectedSuperClusterProducer		=	cms.string('correctedHybridSuperClusters'),
   endcapCorrectedSuperClusterCollection 	=	cms.string(''),
#   endcapCorrectedSuperClusterCollection	=	cms.string('multi5x5EndcapSuperClustersZS'),
   endcapCorrectedSuperClusterProducer		=	cms.string('correctedMulti5x5SuperClustersWithPreshower'),
#   endcapCorrectedSuperClusterProducer         =       cms.string('multi5x5SuperClusters'),
   photonCollection				=	cms.string('photons'),
   electronCollection				= 	cms.string('gsfElectrons'),
   MCParticlesCollection			=   	cms.string('generator'),
   Geant4SimHitsCollection			=	cms.string('g4SimHits'),
   outputFile					=	cms.string('/tmp/hbrun/myOutputFile.root') 
)

#process.HughFilter.ScetThreshold = cms.double(0)
#process.HughFilter.NbXtalThreshold = cms.int32(0)
#process.HughFilter.nEvent = cms.int32(-1)
#process.HughFilter.runNumber = cms.int32(124120)
#process.HughFilter.eventNumbers = cms.vint32(-1)
#process.HughFilter.keepOnlyTwenty = cms.bool(True)


#process.p = cms.Path(process.HughFilter + process.singleEModule)
process.p = cms.Path(process.singleEModule)

