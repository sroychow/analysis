import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  verbosity = cms.untracked.int32(0),
  muonSrc = cms.untracked.InputTag('slimmedMuons'),
  #muonSrc = cms.untracked.InputTag('cleanedMu'),
  #vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  beamSpotCorr = cms.untracked.bool(True),
  muonID = cms.untracked.string('GlobalMuonPromptTight'),
  isMC = cms.untracked.bool(True),
  isSync = cms.untracked.bool(False),
)
