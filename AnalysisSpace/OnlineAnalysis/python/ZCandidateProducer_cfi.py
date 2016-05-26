import FWCore.ParameterSet.Config as cms

zCandidateProducer = cms.EDProducer("ZCandidateProducer",
	verbosity = cms.untracked.int32(0),
	tightIsoElectronFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightEleFSRPairList'),
	tightIsoMuonFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightMuFSRPairList'),
	ZEEColl = cms.untracked.string('ZToeeList'),
	ZMuMuColl = cms.untracked.string('ZTomumuList')
)
