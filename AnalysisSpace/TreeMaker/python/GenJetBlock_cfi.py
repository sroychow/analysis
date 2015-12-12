import FWCore.ParameterSet.Config as cms

genJetBlock = cms.EDProducer('GenJetBlock',
  verbosity = cms.untracked.int32(0),
  genJetSrc = cms.untracked.InputTag('slimmedGenJets')
)
