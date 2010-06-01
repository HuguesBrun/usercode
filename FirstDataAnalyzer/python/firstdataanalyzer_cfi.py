import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("hugues.HughFilter.HughFilter_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'START3X_V26A::All'
process.GlobalTag.globaltag = 'GR_R_35X_V7::All'

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:/tmp/hbrun/theRECOfile.root'
#        '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/133/450/FEEDD762-6C4A-DF11-8A04-003048673E82.root'
'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0185/58E50A5A-6A4E-DF11-BD47-002618943870.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/F6BA6D9E-274E-DF11-875C-002618943885.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/F6394854-274E-DF11-8853-0026189438FF.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/F416A700-284E-DF11-9847-002618943954.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/F2275706-274E-DF11-A46C-002618943896.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/EECED99E-274E-DF11-AF4F-002618943951.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/EE57100C-284E-DF11-A19F-00261894389C.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/EC5ED42A-274E-DF11-ACBC-002618943868.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/E2450806-274E-DF11-85CE-00261894393F.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/E069F104-274E-DF11-9C88-002618943914.root',
#'/store/data/Commissioning10/MinimumBias/RAW-RECO/Apr20Skim_GOODCOLL-v1/0180/D8A54211-274E-DF11-8975-002618943807.root'

)
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#'133450:1-133450:329',
#'133450:332-133450:658'
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.singleEModule = cms.EDAnalyzer('FirstDataAnalyzer',
   readRecHits					=	cms.bool(False),  
   readTTinfos					=	cms.bool(False), 
   readHFrecHits				=	cms.bool(False),
   readMCTruth					=	cms.bool(False),
   readPhotonMCTruth				=	cms.bool(False),
   deltaRMax					= 	cms.double(0.3),
   L1triggerResults 				= 	cms.InputTag('gtDigis'),
   HTLtriggerResults				=	cms.InputTag('TriggerResults','','HLT'),
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
   outputFile					=	cms.string('/tmp/hbrun/myOutputFile_MC.root') 
)

#process.HughFilter.ScetThreshold = cms.double(0)
#process.HughFilter.NbXtalThreshold = cms.int32(0)
#process.HughFilter.nEvent = cms.int32(-1)
#process.HughFilter.runNumber = cms.int32(124120)
#process.HughFilter.eventNumbers = cms.vint32(-1)
#process.HughFilter.keepOnlyTwenty = cms.bool(True)


#process.p = cms.Path(process.HughFilter + process.singleEModule)
process.p = cms.Path(process.singleEModule)

