import FWCore.ParameterSet.Config as cms

process = cms.Process("copy")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("HughFilter_cfi")
process.load('Configuration/EventContent/EventContent_cff')

process.MessageLogger.cerr.threshold = 'INFO'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"file:/sps/cms/sgandurr/CMSSW_4_2_8/src/PhysicsTools/Utilities/scripts/MuMuGammaSelection_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1_AODSIM.root"
)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30)
)


process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands =  process.FEVTEventContent.outputCommands,
                               fileName = cms.untracked.string('MYCOPY.root'),
                               dataset = cms.untracked.PSet(
                                  dataTier = cms.untracked.string('RAW-RECO'),
                                  filterName = cms.untracked.string('Skim_theEvent')),

                               SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring("theHughFilter")
    ))


process.outpath = cms.EndPath(process.out)

process.theHughFilter = cms.Path(process.HughFilter)
#process.p = cms.path(process.HughFilter+process.outpath)


