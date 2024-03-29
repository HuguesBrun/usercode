import FWCore.ParameterSet.Config as cms

process = cms.Process("TotoAna")

# Keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger",
#  cout = cms.untracked.PSet(
#     default = cms.untracked.PSet(
#        limit = cms.untracked.int32(100)
#     ),
#     threshold = cms.untracked.string('INFO')
#   ),
#  destinations = cms.untracked.vstring('cout')
#)

process.options = cms.untracked.PSet(
wantSummary = cms.untracked.bool(True)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Needed for GlobalPositionRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START38_V12::All')
#process.GlobalTag.globaltag = cms.string('GR10_P_V11::All')
#process.GlobalTag.globaltag = cms.string('GR_R_38X_V14::All')

# Global geometry
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')

# Transient Track Builder
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.load("Configuration.StandardSequences.Reconstruction_cff")


# Geometry needed for clustering and calo shapes variables
# process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")

# ES cluster for pi0 discrimination variables
#process.load("RecoEcal.EgammaClusterProducers.preshowerClusterShape_cfi")

# pi0 discrimination variables
#process.load("RecoEcal.EgammaClusterProducers.piZeroDiscriminators_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.maxLuminosityBlocks = cms.untracked.PSet(
#   input = cms.untracked.int32(2)
#)

process.source = cms.Source("PoolSource",

# RECO
fileNames = cms.untracked.vstring('file:/sps/cms/hbrun/dataset_3_8_4_patch3/testRECO/theRECOfile.root')
#,skipEvents=cms.untracked.uint32(26015)
#,firstRun = cms.untracked.uint32(144114),
#,firstLumi   = cms.untracked.uintt32(5),
#,firstEvent = cms.untracked.uint32(135841750)
#eventsToProcess = cms.untracked.VEventRange('144114:135841760-144114:135841771','143192:24655772,-143192:24655772,'),
)

process.totoana = cms.EDAnalyzer("TotoAnalyzer",
   myConfig = cms.PSet(

      # Data type of the PoolSource ( RECO / PAT )
      dataType = cms.untracked.string("RECO"),   # use reco::Objects
      #dataType = cms.untracked.string("PAT"),    # use pat::Objects

      # Verbosite
      #     0 = muet
      #     1 = No evts tous les 10 ou 100 evts
      #     2 = Indique fonctions executees et nb d'objets reconstruits
      #     3 = Liste objets de haut niveau (electrons, muons, photons...)
      #     4 = Liste tous les objets (haut niveau, clusters....)
      #     5 = Debug
      verbosity = cms.untracked.int32(1),

      # name of output root file
      RootFileName = cms.untracked.string('totoOutputMC.root'),

      # DATASET Infos  (will be written in runTree for bookeeping)
      xsection = cms.untracked.double(0.674770994),
      description = cms.untracked.string('Le dataset pourri a Roberto'),

      # What is written to rootuple
      doLHCInfo = cms.untracked.bool(False),
      doL1 = cms.untracked.bool(True),
      doHLT = cms.untracked.bool(True),
      doHLTObject = cms.untracked.bool(True),
      doMC = cms.untracked.bool(True),
      doPDFInfo = cms.untracked.bool(True),
      doSignalMuMuGamma = cms.untracked.bool(False),  # not tested in 2.X.X or 3.X.X
      doSignalTopTop = cms.untracked.bool(False),
      signalGenerator = cms.untracked.string('PYTHIA'),
#     signalGenerator = cms.untracked.string('COMPHEP'),
#     signalGenerator = cms.untracked.string('ALPGEN'),
      #signalGenerator = cms.untracked.string('MADGRAPH'),
      doPhotonConversionMC = cms.untracked.bool(True),

      doPhotonMC = cms.untracked.bool(True),
      doElectronMC = cms.untracked.bool(True),
      doMuonMC = cms.untracked.bool(True),
      doOtherStablePartsMC = cms.untracked.bool(False),
      doJetMC = cms.untracked.bool(True),
      doMETMC = cms.untracked.bool(True),
      doUnstablePartsMC = cms.untracked.bool(True),

      doBeamSpot = cms.untracked.bool(True),
      doPrimaryVertex = cms.untracked.bool(True),
      doZeePrimaryVertex = cms.untracked.bool(True),
      doTrack = cms.untracked.bool(True),
      doJet = cms.untracked.bool(True),
      doMuon = cms.untracked.bool(True),
      doElectron = cms.untracked.bool(True),
      doPhoton = cms.untracked.bool(True),
      doCluster = cms.untracked.bool(True),
      keepClusterizedEcalRecHits = cms.untracked.bool(True),
      keepAllEcalRecHits = cms.untracked.bool(False),
      doMET = cms.untracked.bool(True),
      doBardak = cms.untracked.bool(True),

      doPhotonVertexCorrection = cms.untracked.bool(False),
      doPhotonIsolation = cms.untracked.bool(True),
      doPhotonConversion = cms.untracked.bool(True),
      conversionLikelihoodWeightsFile = cms.untracked.string('RecoEgamma/EgammaTools/data/TMVAnalysis_Likelihood.weights.txt'),

      # Draw MC particle tree
      drawMCTree = cms.untracked.bool(False),
      mcTreePrintP4 = cms.untracked.bool(True),
      mcTreePrintPtEtaPhi = cms.untracked.bool(False),
      mcTreePrintVertex = cms.untracked.bool(False),
      mcTreePrintStatus = cms.untracked.bool(True),
      mcTreePrintIndex = cms.untracked.bool(True),
      mcTreeStatus = cms.untracked.vint32( 1,2,3 ),   # accepted status codes

      # MC particles acceptance cuts
      photonMC_etaMax = cms.double(3.0),
      photonMC_ptMin = cms.double(2.0),
      electronMC_etaMax = cms.double(3.0),
      electronMC_ptMin = cms.double(2.0),
      muonMC_etaMax = cms.double(3.0),
      muonMC_ptMin = cms.double(0.0),
      otherStablePartMC_etaMax = cms.double(3.0),
      otherStablePartMC_ptMin = cms.double(2.0),
      jetMC_etaMax = cms.double(10.0),
      jetMC_ptMin = cms.double(0.0),

      # Photon isolation
      basicClustersIsolation_BarrelBC_type = cms.int32(210),                  # Type of Clusters used for isolation in barrel (see TRootCluster.h for type definition)
      basicClustersIsolation_EndcapBC_type = cms.int32(320),                  # Type of Clusters used for isolation in endcap (see TRootCluster.h for type definition)
      basicClustersIsolation_DRmax = cms.double(0.3),                         # size of the DR cone around photon - Et of BC in this cone are added
      basicClustersIsolation_ClusterEt_threshold = cms.double(0.0),           # Et threshold for BC added in DR cone
      basicClustersDoubleConeIsolation_BarrelBC_type = cms.int32(210),        # Type of Clusters used for isolation in barrel (see TRootCluster.h for type definition)
      basicClustersDoubleConeIsolation_EndcapBC_type = cms.int32(320),        # Type of Clusters used for isolation in endcap (see TRootCluster.h for type definition)
      basicClustersDoubleConeIsolation_DRmin = cms.double(0.05),              # size of the inner DR cone around photon - BC in this cone are rejected
      basicClustersDoubleConeIsolation_DRmax = cms.double(0.3),               # size of the outer DR cone around photon - Et of BC with DRmin < DR < DRmax are added
      basicClustersDoubleConeIsolation_ClusterEt_threshold = cms.double(0.0), # Et threshold for BC added in DR cone
      hcalRecHitIsolation_DRmax = cms.double(0.3),                            # size of the DR cone around photon - Et of HCAL rechits in this cone are added
      hcalRecHitIsolation_HitEt_threshold = cms.double(0.0),                  # Et threshold for HCAL rechits in DR cone
      trackerIsolation_DRmax = cms.double(0.3),                               # size of the DR cone around photon - pt of tracks in this cone are added
      trackerIsolation_pt_threshold = cms.double(0.0),                        # pt threshold for tracks added in DR cone
      trackerIsolation_pixelLayers_threshold = cms.int32(0),                  # minimum number of pixel layers with measurement required for tracks to be added in DR cone isolation

      # Parametrization of the Primary Vertex re-Reconstruction (used for Zee events)
      verbose = cms.untracked.bool(False),
      algorithm = cms.string('AdaptiveVertexFitter'),
      useBeamConstraint = cms.bool(True),
      beamSpotLabel = cms.InputTag("offlineBeamSpot"),
      minNdof  = cms.double(2.0),
      TrackLabel = cms.InputTag("generalTracks"), # label of tracks to be used
                       
      PVSelParameters = cms.PSet(
         maxDistanceToBeam = cms.double(2), ## (in cm) 200 / 500 microns if useBeamConstraint = true / false
         #minVertexFitProb = cms.double(0.01) ## 1% vertex fit probability
      ),

      TkFilterParameters = cms.PSet(
         maxNormalizedChi2 = cms.double(20.0),
         minSiliconLayersWithHits = cms.int32(5), # >= 5
         minPixelLayersWithHits = cms.int32(2),   # >= 2
         maxD0Significance = cms.double(100.0),     # keep most primary tracks
         minPt = cms.double(0.0),                 # better for softish events
         trackQuality = cms.string("any")
      ),

      # clustering
      TkClusParameters = cms.PSet(
         algorithm   = cms.string('gap'),
         TkGapClusParameters = cms.PSet(
            zSeparation = cms.double(0.2) ## 2 mm max separation betw. clusters
         )
      )
      
   ),

   producersNamesRECO = cms.PSet(
      dataType = cms.untracked.string("RECO"),
      allowMissingCollection = cms.untracked.bool(True),
      l1Producer = cms.InputTag("gtDigis"),
      hltProducer = cms.InputTag("TriggerResults","","HLT"),
      hltEvent = cms.InputTag("patTriggerEvent","","HLT"),
      genParticlesProducer = cms.InputTag("genParticles"),
      genJetsProducer = cms.InputTag("antikt5GenJets"),
      genMETsProducer = cms.InputTag("genMetTrue"),
      beamSpotProducer = cms.InputTag("offlineBeamSpot"),
      primaryVertexProducer = cms.InputTag("offlinePrimaryVerticesWithBS"),
      trackProducer = cms.InputTag("generalTracks"),
      jetProducer = cms.VInputTag(
         cms.InputTag("ak5CaloJets"),
         cms.InputTag("kt4PFJets"),
         ),
      muonProducer = cms.VInputTag(cms.InputTag("muons")),
      electronProducer = cms.VInputTag(cms.InputTag("gsfElectrons")),
      photonProducer = cms.VInputTag(cms.InputTag("photons")),
      metProducer = cms.VInputTag(cms.InputTag("met")),
      barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
      endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
      reducedBarrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
      reducedEndcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
      caloTowerCollection = cms.InputTag("towerMaker"),
      hbheRecHitProducer = cms.InputTag("hbhereco"),
      hoRecHitProducer = cms.InputTag("horeco"),
      hfRecHitProducer = cms.InputTag("hfreco"),
      electronProducerForZeeVertex = cms.InputTag("gsfElectrons")
   ),

   producersNamesPAT = cms.PSet(
      dataType = cms.untracked.string("PAT"),
      allowMissingCollection = cms.untracked.bool(True),
      patEncapsulation = cms.untracked.bool(False),
      l1Producer = cms.InputTag("gtDigis"),
      hltProducer = cms.InputTag("TriggerResults","","HLT"),
      hltEvent = cms.InputTag("patTriggerEvent","","PAT"),
      genParticlesProducer = cms.InputTag("genParticles"),
      genJetsProducer = cms.InputTag("antikt5GenJets"),
      genMETsProducer = cms.InputTag("genMetTrue"),
      beamSpotProducer = cms.InputTag("offlineBeamSpot"),
      #primaryVertexProducer = cms.InputTag("offlinePrimaryVerticesWithBS"),
      primaryVertexProducer = cms.InputTag("offlinePrimaryVertices"),
      trackProducer = cms.InputTag("generalTracks"),
      jetProducer = cms.VInputTag(cms.InputTag("cleanPatJets")),
      muonProducer = cms.VInputTag(cms.InputTag("cleanPatMuons")),
      electronProducer = cms.VInputTag(cms.InputTag("cleanPatElectrons")),
      photonProducer = cms.VInputTag(cms.InputTag("cleanPatPhotons")),
      metProducer = cms.VInputTag(cms.InputTag("patMETs")),
      barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
      endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
      reducedBarrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
      reducedEndcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
      caloTowerCollection = cms.InputTag("towerMaker"),
      hbheRecHitProducer = cms.InputTag("hbhereco"),
      hoRecHitProducer = cms.InputTag("horeco"),
      hfRecHitProducer = cms.InputTag("hfreco")
   )
 )


##process.hltHighLevel = cms.EDFilter("HLTHighLevel",
##    HLTPaths = cms.vstring('HLT2PhotonRelaxed'),
##    andOr = cms.bool(True),
##    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
##)

process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
  applyfilter = cms.untracked.bool(True),
  debugOn = cms.untracked.bool(False),
  numtrack = cms.untracked.uint32(10),
  thresh = cms.untracked.double(0.25)
)
import HLTrigger.HLTfilters.triggerResultsFilter_cfi 
process.streamA_datasetPhoton_selector = HLTrigger.HLTfilters.triggerResultsFilter_cfi.triggerResultsFilter.clone()
process.streamA_datasetPhoton_selector.hltResults = cms.InputTag('TriggerResults', '', 'HLT')
process.streamA_datasetPhoton_selector.l1tResults = cms.InputTag('')
process.streamA_datasetPhoton_selector.throw      = cms.bool(False)
process.streamA_datasetPhoton_selector.triggerConditions = cms.vstring('HLT_Photon20_Cleaned_L1R', 
    'HLT_DoublePhoton5_CEP_L1R', 
    'HLT_Photon30_Cleaned_L1R', 
    'HLT_DoublePhoton17_L1R', 
    'HLT_Photon50_NoHE_Cleaned_L1R')



# TotoAna standalone
process.p = cms.Path(process.primaryVertexFilter+process.noscraping+process.streamA_datasetPhoton_selector+process.totoana)
#process.p = cms.Path(process.primaryVertexFilter*process.noscraping*process.streamA_datasetPhoton_selector*process.conversionSequence*process.photonSequence*process.photonIDSequence*process.totoana)

# Photon reReco + TotoAna
#process.load("photonReReco")
#process.p = cms.Path(process.photons*process.photonIDSequence*process.totoana)

# Pi0disc + TotoAna
#process.p = cms.Path(process.preshowerClusterShape*process.piZeroDiscriminators*process.totoana)

# PrimaryVertexFilter + noscraping + hltPhysicsDeclared + TotoAna
#process.p = cms.Path(process.primaryVertexFilter+process.noscraping+process.hltPhysicsDeclared+process.totoana)
