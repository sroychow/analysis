import FWCore.ParameterSet.Config as cms

metBlock = cms.EDProducer('METBlock',
  verbosity = cms.untracked.int32(0),
  metSrc = cms.untracked.InputTag('slimmedMETs'),
  corrmetSrc = cms.untracked.InputTag('slimmedMETsNoHF'), #noHF met temporarily
  #corrmetSrc = cms.untracked.InputTag('patMETsTypeIcorrected'),
  puppimetSrc = cms.untracked.InputTag('slimmedMETsPuppi'),
  mvametSrc = cms.untracked.InputTag('pfMVAMEt')
)
