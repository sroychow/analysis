import FWCore.ParameterSet.Config as cms

eleMVAproducer = cms.EDProducer("ElectronMVAIdProducer",
  verbosity = cms.bool(False),
  electronSrc = cms.InputTag('slimmedElectrons'),
  mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
  mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories")
)

