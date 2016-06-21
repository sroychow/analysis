import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity = cms.untracked.int32(1),
  beamSpotCorr = cms.untracked.bool(True),
  useTrigMode = cms.untracked.bool(False),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  #vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
  electronSrc = cms.untracked.InputTag('slimmedElectrons'),
  useCalibrationMode = cms.untracked.bool(False)  
)

calibelectronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity = cms.untracked.int32(1),
  beamSpotCorr = cms.untracked.bool(True),
  useTrigMode = cms.untracked.bool(False),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  #vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
  electronSrc = cms.untracked.InputTag('calibratedPatElectrons'),
  useCalibrationMode = cms.untracked.bool(True)  
)

