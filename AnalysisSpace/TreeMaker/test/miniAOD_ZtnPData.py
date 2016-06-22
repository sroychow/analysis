import FWCore.ParameterSet.Config as cms
process = cms.Process("TreeMaker")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.infos.threshold = cms.untracked.string("ERROR")
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
  	'/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/02D9C19F-571A-E611-AD8E-02163E013732.root'
  )
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#-----------------------------
# Geometry
#-----------------------------
#process.load("Configuration.StandardSequences.Geometry_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#####mc#####
#process.GlobalTag.globaltag = ''
#####data###
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8'
#process.GlobalTag.globaltag = 'auto:run2_data'
#GR_P_V_56
#
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('ZTnPdata.root')
)
#--------------------------------------------------
# Analysis Tree Specific
#--------------------------------------------------
process.load("AnalysisSpace.TreeMaker.TreeCreator_cfi")
process.load("AnalysisSpace.TreeMaker.TreeWriter_cfi")
process.load("AnalysisSpace.TreeMaker.TreeContentConfig_data_cff")


# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)


########################################################################
#########online analysis sequence#######################
process.load("AnalysisSpace.OnlineAnalysis.TriggerFilter_cfi")
process.load("AnalysisSpace.OnlineAnalysis.ZTnPSelector_cfi")


########################################################

process.p = cms.Path(
  process.selectedPrimaryVertices * 
  process.treeCreator *
  #process.triggerFilterZTnP *
  process.treeContentSequenceZTnPData*
  #process.zTnPanalyzerdata *
  process.treeWriter
)
