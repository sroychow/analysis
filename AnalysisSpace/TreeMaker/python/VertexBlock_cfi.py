import FWCore.ParameterSet.Config as cms

vertexBlock = cms.EDAnalyzer("VertexBlock",
  verbosity = cms.untracked.int32(0),
  #vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices')
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices')
)
