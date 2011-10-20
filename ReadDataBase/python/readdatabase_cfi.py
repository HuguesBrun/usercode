import FWCore.ParameterSet.Config as cms

process = cms.Process("READTHEDB")

# Needed for GlobalPositionRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START42_V11::All')
#process.GlobalTag.globaltag = cms.string('GR_P_V20::All')
#process.GlobalTag.globaltag = cms.string('GR_R_42_V11A::All')
#process.GlobalTag.globaltag = cms.string('GR_P_V18::All')


# Global geometry
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#process.load('Configuration/StandardSequences/GeometryExtended_cff')
#process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'file:/sps/cms/obondu/CMSSW_4_2_3_patch2/src/GluGlu_RECO.root'
)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)


process.readTheDataBase = cms.EDAnalyzer('ReadDataBase'
)

process.p = cms.Path(process.readTheDataBase)
