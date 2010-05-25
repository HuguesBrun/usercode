import FWCore.ParameterSet.Config as cms

process = cms.Process("copy")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("HughFilter_cfi")
process.load('Configuration/EventContent/EventContent_cff')

process.MessageLogger.cerr.threshold = 'INFO'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_Jan29_v8/3dc62ed5c705a8d2a0012a5326b1446e/bscFilter_124275_3.root"
)
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
process.HughFilter.nEvent = cms.int32(1)
process.HughFilter.runNumber = cms.int32(135149)
process.HughFilter.eventNumbers = cms.vint32(125426133)

process.theHughFilter = cms.Path(process.HughFilter)
#process.p = cms.path(process.HughFilter+process.outpath)


