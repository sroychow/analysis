import FWCore.ParameterSet.Config as cms

zzCandidateProducer = cms.EDProducer("ZZCandidateProducer",
	verbosity = cms.untracked.int32(0),
	ZeeSrc = cms.untracked.InputTag('zCandidateProducer','ZToeeList'),
	ZmumuSrc = cms.untracked.InputTag('zCandidateProducer','ZTomumuList'),
	ZZColl = cms.untracked.string('ZZCandList')
)
