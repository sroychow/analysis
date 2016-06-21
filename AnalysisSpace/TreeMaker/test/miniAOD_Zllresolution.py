import FWCore.ParameterSet.Config as cms
process = cms.Process("TreeMaker")
#------------------------
# Message Logger Settings
#------------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.infos.threshold = cms.untracked.string("ERROR")
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/E2490ECF-CBA7-E511-9B19-001E67398458.root',
'file:/afs/cern.ch/user/s/sroychow/public/sinpStorage/sroychow/HZZ4l/syncInput/282C35FB-68A3-E511-A0C4-0CC47A4C8E5E.root',
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
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
#process.GlobalTag.globaltag = 'auto:run2_data'
#GR_P_V_56
#
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('ZllmassResolution.root')
)
#--------------------------------------------------
# Analysis Tree Specific
#--------------------------------------------------
process.load("AnalysisSpace.TreeMaker.TreeCreator_cfi")
process.load("AnalysisSpace.TreeMaker.TreeWriter_cfi")
process.load("AnalysisSpace.TreeMaker.TreeContentConfig_bx25_cff")


# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)

#####Electron Calibrator####


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
	calibratedPatElectrons = cms.PSet(
	initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
        ),
	calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
        ),
)




process.selectedElectrons = cms.EDFilter("PATElectronSelector", 
	src = cms.InputTag("slimmedElectrons"), 
	cut = cms.string("pt > 5 && abs(eta)<2.5") 
)

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedPatElectrons.electrons = "selectedElectrons"
process.calibratedPatElectrons.isMC = True
############################


########################################################################
#########online analysis sequence#######################
process.load("AnalysisSpace.OnlineAnalysis.TriggerFilter_cfi")
process.load("AnalysisSpace.OnlineAnalysis.ZTnPSelector_cfi")


########################################################

process.p = cms.Path(
  process.selectedPrimaryVertices * 
  process.selectedElectrons *
  process.calibratedPatElectrons *
  process.treeCreator *
  process.treeContentSequenceLepCalibZMC*
  process.treeWriter
)
