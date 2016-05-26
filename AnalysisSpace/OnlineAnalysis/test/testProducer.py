import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/mc/RunIIFall15MiniAODv1/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/E2490ECF-CBA7-E511-9B19-001E67398458.root',
#'/store/mc/RunIIFall15MiniAODv1/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/282C35FB-68A3-E511-A0C4-0CC47A4C8E5E.root',
#'/store/mc/RunIIFall15MiniAODv1/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/E2DA5AA7-C5AC-E511-97E0-0CC47A4C8E98.root',

'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/282C35FB-68A3-E511-A0C4-0CC47A4C8E5E.root',
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/E2490ECF-CBA7-E511-9B19-001E67398458.root',
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/E2DA5AA7-C5AC-E511-97E0-0CC47A4C8E98.root',
    )
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'

#--------------------------------------------------
# Analysis Tree Specific
#--------------------------------------------------
#process.load("AnalysisSpace.OnlineAnalysis.HZZFourLElectronSelector_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.HZZFourLMuonSelector_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.FSRBlock_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.IsoLeptonProducer_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.HZZFourLJetSelector_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.ZCandidateProducer_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.ZZCandidateProducer_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.Final4lEventAnalyzer_cfi")
process.load("AnalysisSpace.OnlineAnalysis.TriggerFilter_cfi")
process.load("AnalysisSpace.OnlineAnalysis.OnlineAnalysisMC_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.")

# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)


process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test.root')
)

#process.final4lAnalyzer.computeKD = cms.untracked.bool(True)

#process.p = cms.Path(
#  process.triggerFilter
#  * process.selectedPrimaryVertices 
#  * process.electronSelectorBlock 
#  * process.muonSelectorBlock 
#  * process.fsrBlock
#  * process.isoLeptonProducer
#  * process.zCandidateProducer
#  * process.lepCleanedllooseJetproducer
#  * process.zzCandidateProducer
#  * process.final4lAnalyzer
#)

#process.p = cms.Path(
#	process.triggerFilter
# 	* process.selectedPrimaryVertices
#	* process.onlineAnalysisMCwithKDstandalone
#)

process.p = cms.Path(
	process.triggerFilter
 	* process.selectedPrimaryVertices
	* process.onlineAnalysisMCnoKDstandalone
)


process.out = cms.OutputModule("PoolOutputModule",
     ## Parameters directly for PoolOutputModule
     fileName = cms.untracked.string('test.root'),
)


process.outpath = cms.EndPath(process.out)
