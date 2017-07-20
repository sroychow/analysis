import FWCore.ParameterSet.Config as cms

process = cms.Process("OnlineAnalysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/user/s/sroychow/Exotic/july17/step3miniAOD/Exotic_step3_RAW2DIGI_L1Reco_RECO_EI_PAT_1.root'
    )
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v7'


# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && abs(position.Rho) <= 2"),
  filter = cms.bool(True)                                          
)

###Muon Calibration/Cleaning
# Kalman Muon Calibrations
process.load("AnalysisSpace.OnlineAnalysis.KaMuCaProducer_cfi")
process.calibratedMuons.muonsCollection = cms.InputTag("slimmedMuons")
process.calibratedMuons.isMC = cms.bool(True)
process.calibratedMuons.isSync = cms.bool(True)

# clean muons by segments 
process.gcleanMuons = cms.EDProducer("PATMuonCleanerBySegments",
                                     src = cms.InputTag("calibratedMuons"),
                                     preselection = cms.string("track.isNonnull"),
                                     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                     fractionOfSharedSegments = cms.double(0.499),
                                     )

###Electron Calibration Cleaning, MV ID#####################
# Electron Smear/Scale
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         cut = cms.string("pt > 5 && abs(eta)<2.5")
                                         )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(123456),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
process.calibratedPatElectrons.electrons = cms.InputTag('selectedElectrons')
process.calibratedPatElectrons.isMC = cms.bool(True)
process.calibratedPatElectrons.isSynchronization = cms.bool(False)
process.calibratedPatElectrons.correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_Golden22June_approval")

##
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
fileFormat = 'miniAOD'

if fileFormat == 'AOD' :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)


# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("calibratedPatElectrons")


###Jet part####
## Sequence : 1. JEC 2. QGtagger 3. Producer with QG+JER from JEC applied jets
##############
#############################################################################
#ReApply JEC on Data
#############################################################################
#Remember to change the jet collection name as 'patJetsReapplyJEC'
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors

process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    ##levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], ##for DATA
    levels = ['L1FastJet', 'L2Relative', 'L3Absolute'], ##forMC
    payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets

process.patJetsReapplyJEC = updatedPatJets.clone(
      jetSource = cms.InputTag("slimmedJets"),
      jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
      )

##############################################################################
##QG tagger
##############################################################################
process.load("CondCore.CondDB.CondDB_cfi")
qgDatabaseVersion = 'cmssw8020_v2'
import os
QGddatabasefile = os.environ.get('CMSSW_BASE')+"/src/AnalysisSpace/OnlineAnalysis/data/QGL_"+qgDatabaseVersion+".db"
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('sqlite_file:' + QGddatabasefile)
)
process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

###########JER + QG producer#######
process.load("AnalysisSpace.OnlineAnalysis.JetQGtaggProducer_cfi")
process.skimmedJetswqg.jetSrc = cms.InputTag('patJetsReapplyJEC')
#########################MET###################################
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
                           isData=True,#(or False),
                           )
"""
# If you would like to re-cluster and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
                           isData=True, # or False),
                           pfCandColl=cms.InputTag("packedPFCandidates"),
                           recoMetFromPFCs=True,
                           )

# If you would like to re-cluster both jets and met and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process,
                           isData=True,#(or False),
                           pfCandColl=cms.InputTag("packedPFCandidates"),
                           recoMetFromPFCs=True,
                           CHS = True, #This is an important step and determines what type of jets to be reclustered
                           reclusterJets = True
                           )
"""
######################MET###############################################
##############################################################################
#---------------------------------------------------------------------------
# Trigger filter and EleMVaId producer from Online
#---------------------------------------------------------------------------
process.load("AnalysisSpace.OnlineAnalysis.TriggerFilter_cfi")
process.load("AnalysisSpace.OnlineAnalysis.ElectronMVAIdProducer_cfi")
process.eleMVAproducer.electronSrc = cms.InputTag('calibratedPatElectrons')

####TreeMaker Part####
process.load("AnalysisSpace.TreeMaker.TreeCreator_cfi")
process.load("AnalysisSpace.TreeMaker.TreeWriter_cfi")
process.load("AnalysisSpace.TreeMaker.TreeContentConfig_bx25_cff")
##TreeMaker will take calibrated electrons and muons
process.muonBlock.muonSrc = cms.InputTag('gcleanMuons') 
process.electronBlock.electronSrc = cms.InputTag('eleMVAproducer')




process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('test.root')
                                   )


process.p = cms.Path(process.triggerFilter
 	             * process.selectedPrimaryVertices
                     * process.calibratedMuons
                     * process.gcleanMuons
                     * process.regressionApplication
                     * process.selectedElectrons
                     * process.calibratedPatElectrons
                     * process.electronMVAValueMapProducer
                     * process.eleMVAproducer
                     * process.patJetCorrFactorsReapplyJEC
                     * process.patJetsReapplyJEC
                     * process.QGTagger
                     * process.skimmedJetswqg
                     * process.fullPatMetSequence
                     * process.treeCreator
                     * process.treeContentSequence
                     * process.treeWriter
                     )
