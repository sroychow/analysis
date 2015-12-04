import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/mc/RunIISpring15MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/3E964C5D-1D6E-E511-8B9A-0050560207C5.root',
    #'/store/mc/RunIISpring15MiniAODv2/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/D8CA6B54-056F-E511-BB1A-02163E014CE3.root',
    #'/store/mc/RunIISpring15MiniAODv2/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/D22BEE88-C26D-E511-B330-002590A81EF0.root',
    '/store/mc/RunIISpring15MiniAODv2/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/80000/84F62DD7-1475-E511-9F59-009C02AB98A6.root'
    )
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'


#--------------------------------------------------
# Analysis Tree Specific
#--------------------------------------------------
process.load("AnalysisSpace.TreeMaker.TreeCreator_cfi")
process.load("AnalysisSpace.TreeMaker.TreeWriter_cfi")
#process.load("AnalysisSpace.TreeMaker.TreeContentConfig_data_cff")
process.load("AnalysisSpace.TreeMaker.TreeContentConfig_bx25_cff")

process.calibratedElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    electrons = cms.InputTag("slimmedElectrons"),
    grbForestName = cms.string("gedelectron_p4combination_25ns"),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(True)
)


# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
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

process.egmGsfElectronIDsSlimmed = process.egmGsfElectronIDs.clone(physicsObjectSrc = cms.InputTag("slimmedElectrons"))
process.electronMVAValueMapProducerSlimmed = process.electronMVAValueMapProducer.clone(
                                           srcMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"))
process.electronRegressionValueMapProducerSlimmed = process.electronRegressionValueMapProducer.clone(
                                             srcMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"))
process.egmGsfElectronIDSequenceSlimmed = cms.Sequence( process.electronMVAValueMapProducerSlimmed +
                                                       process.egmGsfElectronIDsSlimmed +
                                                       process.electronRegressionValueMapProducerSlimmed)

process.eleMVASequence = cms.Sequence( process.egmGsfElectronIDSequenceSlimmed + process.egmGsfElectronIDSequence )


process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('calibratedElectrons')

########################################################################
process.path = cms.Path(
    process.calibratedElectrons*
    process.eleMVASequence*
    process.selectedPrimaryVertices*
    process.treeCreator*
    process.treeContentSequence*
    process.treeWriter
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('HZZ4lminiAOD_CATSYNC2ndDec.root'))

