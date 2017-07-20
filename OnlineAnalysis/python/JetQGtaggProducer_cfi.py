import FWCore.ParameterSet.Config as cms
skimmedJetswqg = cms.EDProducer("JetQGproducer",
  verbosity = cms.bool(False),
  isMC = cms.bool(True),
  jetPtcut = cms.double(30.),
  jetEtacut = cms.double(4.7),
  muRho = cms.InputTag("fixedGridRhoFastjetAll"),
  jetSrc = cms.InputTag('slimmedJets') 
)
