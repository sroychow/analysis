import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDProducer("ElectronBlock",
  verbosity = cms.untracked.int32(1),
  beamSpotCorr = cms.untracked.bool(True),
  useTrigMode = cms.untracked.bool(False),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  #vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
  electronSrc = cms.untracked.InputTag('slimmedElectrons'),
  eleMediumIdMap = cms.InputTag("egmGsfElectronIDsSlimmed:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
  eleTightIdMap = cms.InputTag("egmGsfElectronIDsSlimmed:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),  
  mvaValuesMap     = cms.InputTag("electronMVAValueMapProducerSlimmed:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
  mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducerSlimmed:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
  objectbranchName = cms.string("Electron")
)

