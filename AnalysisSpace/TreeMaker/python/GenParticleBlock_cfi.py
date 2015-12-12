import FWCore.ParameterSet.Config as cms

genParticleBlock = cms.EDProducer("GenParticleBlock",
  verbosity = cms.untracked.int32(0),
  genParticleSrc = cms.untracked.InputTag('prunedGenParticles')
)
