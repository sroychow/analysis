import FWCore.ParameterSet.Config as cms

from AnalysisSpace.OnlineAnalysis.HZZFourLElectronSelector_cfi import *
from AnalysisSpace.OnlineAnalysis.HZZFourLMuonSelector_cfi import *
from AnalysisSpace.OnlineAnalysis.FSRBlock_cfi import *
from AnalysisSpace.OnlineAnalysis.IsoLeptonProducer_cfi import *
from AnalysisSpace.OnlineAnalysis.HZZFourLJetSelector_cfi import * 
from AnalysisSpace.OnlineAnalysis.ZCandidateProducer_cfi import *
from AnalysisSpace.OnlineAnalysis.ZZCandidateProducer_cfi import *
from AnalysisSpace.OnlineAnalysis.Final4lEventAnalyzerMC_cfi import *

onlineSelector = cms.Sequence(
    muonSelectorBlock 
  * electronSelectorBlock
  * fsrBlock
  * isoLeptonProducer
  * zCandidateProducer
  * lepCleanedllooseJetproducer
  * zzCandidateProducer
)
####### MC with KD & standalone###########
onlineAnalysisMCwithKDstandalone = cms.Sequence(onlineSelector * final4lAnalyzerwithKDstandalone)

####### MC without KD & standalone###########
onlineAnalysisMCnoKDstandalone = cms.Sequence(onlineSelector * final4lAnalyzernoKDstandalone)

####### MC with KD & savetree###########
onlineAnalysisMCwithKDsaveEvent = cms.Sequence(onlineSelector * final4lAnalyzerwithKDsaveTree)
