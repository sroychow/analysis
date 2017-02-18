import FWCore.ParameterSet.Config as cms

process = cms.Process("OnlineAnalysis")

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
    '/store/mc/RunIISummer16MiniAODv2/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22F32262-3FC5-E611-B373-D4AE526DEDB7.root',
    '/store/mc/RunIISummer16MiniAODv2/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/387FA719-E6CC-E611-A1F0-FA163E7D6032.root',
    '/store/mc/RunIISummer16MiniAODv2/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/50DE4DA2-1EC1-E611-9A3C-002590E7E010.root'
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

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
                                        # input collections
                                        electrons = cms.InputTag('selectedElectrons'),
                                        gbrForestName = cms.string("gedelectron_p4combination_25ns"),
                                        isMC = cms.bool(True),
                                        isSynchronization = cms.bool(True),
                                        correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/80X_Golden22June_approval")
                                        )
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
my_id_modules = [
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_V1_cff',
                ]

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

##############################################################################
#---------------------------------------------------------------------------
# Analysis Tree Specific
#---------------------------------------------------------------------------
process.load("AnalysisSpace.OnlineAnalysis.TriggerFilter_cfi")
process.load("AnalysisSpace.OnlineAnalysis.FourLeptonSelector_cfi")
process.load("AnalysisSpace.OnlineAnalysis.ElectronMVAIdProducer_cfi")
#process.load("AnalysisSpace.OnlineAnalysis.")

process.eleMVAproducer.electronSrc = cms.InputTag('calibratedPatElectrons')
#process.eleMVAproducer.electronSrc = cms.InputTag('calibratedPatElectrons')
process.fourleptonSelector.muonSrc = cms.InputTag('gcleanMuons') 
process.fourleptonSelector.electronSrc = cms.InputTag('eleMVAproducer')


process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test.root')
)


process.p = cms.Path(process.triggerFilter
 	             * process.selectedPrimaryVertices
                     * process.calibratedMuons
                     * process.gcleanMuons
                     * process.selectedElectrons
                     * process.calibratedPatElectrons
                     * process.electronMVAValueMapProducer
                     * process.eleMVAproducer
                     * process.patJetCorrFactorsReapplyJEC
                     * process.patJetsReapplyJEC
                     * process.QGTagger
                     * process.skimmedJetswqg
                     * process.fourleptonSelector
)
