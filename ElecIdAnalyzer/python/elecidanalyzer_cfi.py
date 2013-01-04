import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V6::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#
# the MC global Tag : START52_V9
# the DATA gloval Tag : GR_P_V32

# input
#

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
	#'file:/sps/cms/hbrun/CMSSW_5_3_2_myCode/src/files/RelValZEE.root'
	#'file:/sps/cms/hbrun/CMSSW_5_3_2_myCode/src/files/DYJetsToLLM50.root'
    #'file:/sps/cms/hbrun/CMSSW_5_3_2_myCode/src/files/mumuFiles/mumugamma_0.root'
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_1_1_5Q6.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_2_1_gpy.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_3_1_lfX.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_4_1_eOA.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_5_1_Boi.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_6_1_KgO.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_7_1_YG7.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_8_1_GeK.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_9_1_5fE.root',
                                      'file:/sps/cms/hbrun/CMSSW_5_2_5_recup/src/files/MYCOPY_10_1_bd2.root'
    ),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
# all the analyses
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

#compute the pf iso
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

process.theDiElecFilter = cms.EDFilter('DiElecFilter',
     electronsInputTag       = cms.InputTag("gsfElectrons")                                       
)
                                    

process.theEleIdAnalyzer = cms.EDAnalyzer('ElecIdAnalyzer',
    isMC                    = cms.bool(True),
	doMuons					= cms.bool(False),
	doPhotons				= cms.bool(False),
	savePF					= cms.bool(True),
	saveConversions			= cms.bool(False),
    doMuMuGammaMC           = cms.bool(False),
    electronsInputTag       = cms.InputTag("gsfElectrons"),
    conversionsInputTag     = cms.InputTag("allConversions"),
    beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
    rhoIsoInputTag          = cms.InputTag("kt6PFJetsForIsolation", "rho"),
    primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
	muonProducer = cms.VInputTag(cms.InputTag("muons")),
	isoValInputTags         = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
											cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
											cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'),
											cms.InputTag('elPFIsoValueChargedAll03PFIdPFIso'),
											cms.InputTag('elPFIsoValuePU03PFIdPFIso'),
											cms.InputTag('elPFIsoValueCharged04PFIdPFIso'),
											cms.InputTag('elPFIsoValueGamma04PFIdPFIso'),
											cms.InputTag('elPFIsoValueNeutral04PFIdPFIso'),	
											cms.InputTag('elPFIsoValueChargedAll04PFIdPFIso'),
											cms.InputTag('elPFIsoValuePU04PFIdPFIso')),
    TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),
    HLTTriggerSummaryAOD    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
	photonCollection		= cms.string("photons"),
    outputFile		        = cms.string("elecIDtree.root"),
    deltaRsavePF            = cms.double(0.6),
    printDebug              = cms.bool(True)
)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)	
                                           )


process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.electronsMCmatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                                     src = cms.InputTag("gsfElectrons"),
                                     matched = cms.InputTag("genParticles"),
                                     distMin = cms.double(0.15),
                                     matchPDGId = cms.vint32( 11 ) # electrons
                                     )

#process.load("hugues.ElecIdAnalyzer.ZmmgSkim_cfi")
#process.myZmmgSkimSeq = cms.Sequence(process.ZmmgSkimSeq)





process.p = cms.Path(process.primaryVertexFilter * process.noscraping * process.kt6PFJetsForIsolation * process.pfiso * process.theEleIdAnalyzer)
