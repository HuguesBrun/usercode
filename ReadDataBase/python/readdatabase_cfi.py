import FWCore.ParameterSet.Config as cms

process = cms.Process("READTHEDB")

# Needed for GlobalPositionRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('FT_R_39X_V4A::All')

# Global geometry
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load('Configuration/StandardSequences/MagneticField_38T_cff')
#process.load('Configuration/StandardSequences/GeometryExtended_cff')
#process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'file:/sps/cms/hbrun/dataset_3_9_7/theRECOfile.root'
)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)


process.readTheDataBase = cms.EDAnalyzer('ReadDataBase'
)

process.p = cms.Path(process.readTheDataBase)
