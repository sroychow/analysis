import FWCore.ParameterSet.Config as cms

###############MC+withKD+Standalone##########################

final4lAnalyzerwithKDstandalone = cms.EDAnalyzer('Final4lEventAnalyzer',
	ZZcandSrc = cms.untracked.InputTag('zzCandidateProducer','ZZCandList'),
	lepCleanedJetSrc = cms.untracked.InputTag('lepCleanedllooseJetproducer','lepCleanedllooseJetList'),
	tightIsoElectronFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightEleFSRPairList'),
	tightIsoMuonFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightMuFSRPairList'),
        syncFilename = cms.untracked.string('syncDump_withKD.txt'),
        runStandalone = cms.untracked.bool(True),
        computeKD = cms.untracked.bool(True),
        isData = cms.untracked.bool(False)
)

###############MC+ withoutKD + Standalone##########################

final4lAnalyzernoKDstandalone = cms.EDAnalyzer('Final4lEventAnalyzer',
	ZZcandSrc = cms.untracked.InputTag('zzCandidateProducer','ZZCandList'),
	lepCleanedJetSrc = cms.untracked.InputTag('lepCleanedllooseJetproducer','lepCleanedllooseJetList'),
	tightIsoElectronFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightEleFSRPairList'),
	tightIsoMuonFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightMuFSRPairList'),
        syncFilename = cms.untracked.string('syncDump_noKD.txt'),
        runStandalone = cms.untracked.bool(True),
        computeKD = cms.untracked.bool(False),
        isData = cms.untracked.bool(False)
)

###############MC+withKD+saveinTree##########################

final4lAnalyzerwithKDsaveTree = cms.EDAnalyzer('Final4lEventAnalyzer',
	ZZcandSrc = cms.untracked.InputTag('zzCandidateProducer','ZZCandList'),
	lepCleanedJetSrc = cms.untracked.InputTag('lepCleanedllooseJetproducer','lepCleanedllooseJetList'),
	tightIsoElectronFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightEleFSRPairList'),
	tightIsoMuonFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightMuFSRPairList'),
        syncFilename = cms.untracked.string('syncDump_withKD.txt'),
        runStandalone = cms.untracked.bool(False),
        computeKD = cms.untracked.bool(True),
        isData = cms.untracked.bool(False)
)

