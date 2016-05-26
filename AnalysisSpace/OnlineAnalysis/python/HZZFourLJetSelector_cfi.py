import FWCore.ParameterSet.Config as cms

lepCleanedllooseJetproducer = cms.EDProducer("HZZFourLJetSelector",
  verbosity = cms.untracked.int32(0),
  jetSrc = cms.untracked.InputTag('slimmedJets'),
  tightIsoElectronFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightEleFSRPairList'),
  tightIsoMuonFSRSrc = cms.untracked.InputTag('isoLeptonProducer','tightMuFSRPairList'),
  looseJetColl = cms.untracked.string("lepCleanedllooseJetList"),
  jecUncertainty = cms.string('CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt'),
  applyResidualJEC = cms.bool(False),
  residualJEC = cms.string('CondFormats/JetMETObjects/data/Spring10_L2L3Residual_AK5PF.txt')
)
