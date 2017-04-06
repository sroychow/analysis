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
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_1.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_2.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_3.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_4.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_5.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_6.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_7.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_8.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_9.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_10.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_11.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_12.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_13.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_14.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_15.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_16.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_17.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_18.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_19.root',
        'file:/afs/cern.ch/work/s/subir/public/data/Exotic/pph1_170_h2_70/step4/Exotic_MINIAODSIM_20.root'
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
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
#
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Exotic_pph1_170_h2_70_1.root')
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

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs we want to produce
my_id_modules = [
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_V1_cff',
                ]
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.p = cms.Path(
  process.electronMVAValueMapProducer*
  process.selectedPrimaryVertices* 
  process.treeCreator*
  process.treeContentSequence*
  process.treeWriter
)
