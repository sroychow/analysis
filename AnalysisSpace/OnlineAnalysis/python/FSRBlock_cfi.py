import FWCore.ParameterSet.Config as cms

fsrBlock = cms.EDProducer("FSRBlock",
	verbosity = cms.untracked.int32(0),
	pfCands = cms.untracked.InputTag('packedPFCandidates'),
 	looseElectronSrc = cms.untracked.InputTag('electronSelectorBlock','looseSIPElectronVector'),
	looseMuonSrc = cms.untracked.InputTag('muonSelectorBlock','looseSIPMuonVector'),
        selectedFSRColl =  cms.untracked.string('selectedFSRVector'),
	looseSIPEleFSRColl = cms.untracked.string('looseEleFSRPairList'),	
	looseSIPMuFSRColl = cms.untracked.string('looseMuFSRPairList')
)
