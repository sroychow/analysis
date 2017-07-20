import FWCore.ParameterSet.Config as cms

calibratedMuons = cms.EDProducer("KaMuCaProducer",
  verbosity = cms.bool(False),
  isMC = cms.bool(True),
  isSync = cms.bool(False),
  muonSrc = cms.InputTag('slimmedMuons') 
)
