import FWCore.ParameterSet.Config as cms
process = cms.Process("TreeMaker")
#------------------------
# Message Logger Settings
#------------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.infos.threshold = cms.untracked.string("ERROR")
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/mc/RunIISpring15MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/3E964C5D-1D6E-E511-8B9A-0050560207C5.root',
    #'/store/mc/RunIISpring15MiniAODv2/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/D8CA6B54-056F-E511-BB1A-02163E014CE3.root',
    #'/store/mc/RunIISpring15MiniAODv2/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/D22BEE88-C26D-E511-B330-002590A81EF0.root',
    #'/store/mc/RunIISpring15MiniAODv2/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/80000/84F62DD7-1475-E511-9F59-009C02AB98A6.root'
  )
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-----------------------------
# Geometry
#-----------------------------
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
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

# Sequence for rhoNeutral
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi import pfAllNeutralHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import pfAllPhotons
pfNeutralCandPdgIds = []
pfNeutralCandPdgIds.extend(pfAllNeutralHadrons.pdgId.value())
pfNeutralCandPdgIds.extend(pfAllPhotons.pdgId.value())

process.pfNeutralCands = cms.EDFilter("PdgIdPFCandidateSelector",
  src = cms.InputTag('particleFlow'),
  pdgId = cms.vint32(pfNeutralCandPdgIds)
)
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFNeutralJetsForVtxMultReweighting = kt4PFJets.clone(
  src = cms.InputTag('pfNeutralCands'),
  rParam = cms.double(0.6),
  doRhoFastjet = cms.bool(True),
  Rho_EtaMax = cms.double(2.5)
)
process.kt6PFJets = kt4PFJets.clone(
  rParam = cms.double(0.6),
  doRhoFastjet = cms.bool(True),
  Rho_EtaMax = cms.double(2.5)
)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
   src = cms.InputTag("calibratedMuons"),
   preselection = cms.string("track.isNonnull"),
   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
   fractionOfSharedSegments = cms.double(0.499)
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
##############Electron Calibration#####################################
#random number generation
import random
process.load('IOMC.RandomEngine.IOMC_cff')
random.seed()
process.RandomNumberGeneratorService.generator.initialSeed  = random.randrange(1,10e07)

process.selectedSlimmedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt > 7."),
)
process.calibratedElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    electrons = cms.InputTag("selectedSlimmedElectrons"),
    grbForestName = cms.string("gedelectron_p4combination_25ns"),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(False),
)
#######################################################################

process.p = cms.Path(
  ( process.selectedSlimmedElectrons + process.calibratedElectrons ) *
  process.egmGsfElectronIDSequence*
  process.selectedPrimaryVertices * 
  #process.cleanedMu *
  #process.kt6PFJets* 
  #process.pfNeutralCands *
  #process.kt6PFNeutralJetsForVtxMultReweighting*
  process.treeCreator*
  process.treeContentSequence*
  process.treeWriter
)

# List File names here
#---------------------------------------
#process.PoolSource.fileNames = [
#  '/store/mc/Spring14miniaod/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0881ABEB-2709-E411-9E42-00145EDD7581.root'
#]
#process.PoolSource.fileNames = ['']
