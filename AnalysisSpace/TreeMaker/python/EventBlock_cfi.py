import FWCore.ParameterSet.Config as cms

eventBlock = cms.EDAnalyzer("EventBlock",
  verbosity = cms.untracked.int32(0),
  l1InputTag = cms.untracked.InputTag('gtDigis'),
  vertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  vertexMinimumNDOF = cms.untracked.uint32(4),
  vertexMaxAbsZ = cms.untracked.double(24.),
  vertexMaxd0 = cms.untracked.double(2.),
  hpTrackThreshold = cms.untracked.double(0.25),
  puSummaryInputTag = cms.untracked.InputTag('slimmedAddPileupInfo'),
  selectedVtxInputTag = cms.untracked.InputTag('selectedPrimaryVertices'),
  rhoInputTag = cms.untracked.InputTag('kt6PFJets','rho'),                         
  rhoNeutralInputTag = cms.untracked.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho')                         
)
