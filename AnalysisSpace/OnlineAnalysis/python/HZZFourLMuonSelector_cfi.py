import FWCore.ParameterSet.Config as cms

muonSelectorBlock = cms.EDProducer("HZZFourLMuonSelector",
  verbosity = cms.untracked.int32(0),
  muonSrc = cms.untracked.InputTag('slimmedMuons'),
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
  looseMucoll = cms.untracked.string('looseMuonVector'),
  looseSIPMucoll = cms.untracked.string('looseSIPMuonVector'),
  tightMucoll = cms.untracked.string('tightMuonVector')
)
