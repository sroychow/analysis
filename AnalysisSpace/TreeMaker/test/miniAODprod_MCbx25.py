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
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/282C35FB-68A3-E511-A0C4-0CC47A4C8E5E.root',
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/E2490ECF-CBA7-E511-9B19-001E67398458.root',
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/E2DA5AA7-C5AC-E511-97E0-0CC47A4C8E98.root',
  )
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
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
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'
#process.GlobalTag.globaltag = 'auto:run2_data'
#GR_P_V_56
#
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('HZZ4lminiAOD_CATSYNC23rdOct_new.root')
)
#--------------------------------------------------
# Analysis Tree Specific
#--------------------------------------------------
process.load("AnalysisSpace.TreeMaker.TreeCreator_cfi")
process.load("AnalysisSpace.TreeMaker.TreeWriter_cfi")
#process.load("AnalysisSpace.TreeMaker.TreeContentConfig_bx50_cff")
process.load("AnalysisSpace.TreeMaker.TreeContentConfig_bx25_cff")


# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)


#################################################################
#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
########################################################################
#########online analysis sequence#######################
process.load("AnalysisSpace.OnlineAnalysis.TriggerFilter_cfi")
process.load("AnalysisSpace.OnlineAnalysis.OnlineAnalysisMC_cfi")
#process.onlineAnalysis = cms.Sequence(process.triggerFilter*process.onlineAnalysisMCwithKDstandalone)
process.onlineAnalysis = cms.Sequence(process.triggerFilter*process.onlineAnalysisMCwithKDsaveEvent)

########################################################



process.p = cms.Path(
  process.selectedPrimaryVertices * 
  process.treeCreator*
  (process.onlineAnalysis + process.treeContentSequence)*
  process.treeWriter
)

# List File names here
#---------------------------------------
#process.PoolSource.fileNames = [
#  '/store/mc/Spring14miniaod/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0881ABEB-2709-E411-9E42-00145EDD7581.root'
#]
#process.PoolSource.fileNames = ['']
