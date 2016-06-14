import FWCore.ParameterSet.Config as cms

calibMuonprodMCnosync = cms.EDProducer("KalmanMuonCalibrationsProducer",
  muonSrc = cms.InputTag('slimmedMuons'),
  isMC = cms.bool(True),
  isSync = cms.bool(False),
  calibratedMucoll = cms.string('kalmanCalibMuonVector')
)

calibMuonprodMCsync = cms.EDProducer("KalmanMuonCalibrationsProducer",
  muonSrc = cms.InputTag('slimmedMuons'),
  isMC = cms.bool(True),
  isSync = cms.bool(True),
  calibratedMucoll = cms.string('kalmanCalibMuonVector')
)

calibMuonprodData = cms.EDProducer("KalmanMuonCalibrationsProducer",
  muonSrc = cms.InputTag('slimmedMuons'),
  isMC = cms.bool(False),
  isSync = cms.bool(False),
  calibratedMucoll = cms.string('kalmanCalibMuonVector')
)
