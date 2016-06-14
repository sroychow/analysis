import FWCore.ParameterSet.Config as cms

from AnalysisSpace.OnlineAnalysis.HZZFourLElectronSelector_cfi import *
from AnalysisSpace.OnlineAnalysis.HZZFourLMuonSelector_cfi import *
from AnalysisSpace.OnlineAnalysis.FSRBlock_cfi import *

ztnpwithFSRmc = cms.EDAnalyzer("ZTnpAnalyzer",
	verbosity = cms.untracked.int32(0),
	fsrCol = cms.untracked.InputTag('fsrBlock','selectedFSRVector'), 
        vertexSrc = cms.untracked.InputTag('selectedPrimaryVertices'),
	looseSIPElectronFSRSrc = cms.untracked.InputTag('fsrBlock','looseEleFSRPairList'),
	looseSIPMuonFSRSrc = cms.untracked.InputTag('fsrBlock','looseMuFSRPairList'),
	fixedGridRhoFastjetAllTag = cms.untracked.InputTag('fixedGridRhoFastjetAll'),
        isData=cms.untracked.bool(False)
)

zTnPanalyzermc = cms.Sequence(
  muonSelectorBlock 
  * electronSelectorBlock 
  * fsrBlock
  * ztnpwithFSRmc
)

