import FWCore.ParameterSet.Config as cms

fourleptonSelector = cms.EDAnalyzer('FourLeptonSelector',
    verbosity = cms.bool(False),
    isMC = cms.bool(True),
    syncFilename = cms.untracked.string("sync_80X_feb2017.txt"),
    fixedGridRhoFastjetAllTag = cms.untracked.InputTag('fixedGridRhoFastjetAll'),
    vertexSrc = cms.InputTag('selectedPrimaryVertices'),
    muonSrc = cms.InputTag('slimmedMuons'),
    electronSrc = cms.InputTag('slimmedElectrons'),
    jetSrc = cms.InputTag('slimmedJets'),
    pfcandSrc = cms.InputTag('packedPFCandidates'),
)
