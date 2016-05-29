import FWCore.ParameterSet.Config as cms

final4lAnalyzer = cms.EDAnalyzer('Final4lEventAnalyzer',
	ZZcandSrc = cms.untracked.InputTag('zzCandidateProducer','ZZCandList'),
	lepCleanedJetSrc = cms.untracked.InputTag('lepCleanedllooseJetproducer','lepCleanedllooseJetList'),
	tightIsoElectronFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightEleFSRPairList'),
	tightIsoMuonFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightMuFSRPairList'),
        syncFilename = cms.untracked.string('syncDump.txt'),
        runStandalone = cms.untracked.bool(True),
        computeKD = cms.untracked.bool(False),
        isData = cms.untracked.bool(False)
)
