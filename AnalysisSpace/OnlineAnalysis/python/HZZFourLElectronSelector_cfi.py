import FWCore.ParameterSet.Config as cms

electronSelectorBlock = cms.EDProducer("HZZFourLElectronSelector",
  verbosity = cms.untracked.int32(1),
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
  electronSrc = cms.untracked.InputTag('slimmedElectrons'),
  looseElecoll = cms.untracked.string('looseElectronVector'),
  looseSIPElecoll = cms.untracked.string('looseSIPElectronVector'),
  tightElecoll = cms.untracked.string('tightElectronVector')
)

