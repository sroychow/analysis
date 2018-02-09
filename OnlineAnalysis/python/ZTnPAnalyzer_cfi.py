import FWCore.ParameterSet.Config as cms

zTnPAnalyzer = cms.EDAnalyzer('ZTnpAnalyzer',
    verbosity = cms.bool(False),
    isMC = cms.bool(True),
    fixedGridRhoFastjetAllTag = cms.untracked.InputTag('fixedGridRhoFastjetAll'),
    vertexSrc = cms.InputTag('selectedPrimaryVertices'),
    muonSrc = cms.InputTag('slimmedMuons'),
    electronSrc = cms.InputTag('slimmedElectrons'),
    pfcandSrc = cms.InputTag('packedPFCandidates'),
)
