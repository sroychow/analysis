import FWCore.ParameterSet.Config as cms

genMETBlock = cms.EDProducer('GenMETBlock',
  verbosity = cms.untracked.int32(0),
  genMETSrc = cms.untracked.InputTag('genMetTrue')
)
