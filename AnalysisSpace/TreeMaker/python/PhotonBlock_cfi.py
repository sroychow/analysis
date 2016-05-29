import FWCore.ParameterSet.Config as cms

photonBlock = cms.EDAnalyzer("PhotonBlock",
  verbosity = cms.untracked.int32(1),
  photonSrc = cms.untracked.InputTag('slimmedPhotons')
)
