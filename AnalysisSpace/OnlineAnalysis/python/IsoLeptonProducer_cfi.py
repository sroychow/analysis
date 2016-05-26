import FWCore.ParameterSet.Config as cms

isoLeptonProducer = cms.EDProducer("IsoLeptonProducer",
	verbosity = cms.untracked.int32(0),
	fsrCol = cms.untracked.InputTag('fsrBlock','selectedFSRVector'), 
	looseSIPElectronFSRSrc = cms.untracked.InputTag('fsrBlock','looseEleFSRPairList'),
	looseSIPMuonFSRSrc = cms.untracked.InputTag('fsrBlock','looseMuFSRPairList'),
	fixedGridRhoFastjetAllTag = cms.untracked.InputTag('fixedGridRhoFastjetAll'),
	tightSIPEleFSRColl = cms.untracked.string('tightEleFSRPairList'),
	tightSIPMuFSRColl = cms.untracked.string('tightMuFSRPairList')
)

