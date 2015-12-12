import FWCore.ParameterSet.Config as cms

packedPFCandidateBlock = cms.EDProducer("PackedPFCandidateBlock",
  verbosity = cms.untracked.int32(0),
  pfCands   = cms.untracked.InputTag('packedPFCandidates'),
  pdgTosave = cms.vint32( 22 )  
)
