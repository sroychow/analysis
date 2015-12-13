#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/METBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

METBlock::METBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc", edm::InputTag("patMETs"))),
  corrMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("corrmetSrc", edm::InputTag("patMETsTypeIcorrected"))),
  puppiMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("puppimetSrc", edm::InputTag("patMETsTypeIcorrected"))),
//  mvaMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvametSrc", edm::InputTag("patPFMetByMVA"))),
  mvaMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvametSrc", edm::InputTag("pfMVAMEt"))),
  pfMETToken_(consumes<pat::METCollection>(pfMETTag_)),
  corrMETToken_(consumes<pat::METCollection>(corrMETTag_)),
  puppiMETToken_(consumes<pat::METCollection>(puppiMETTag_)),
//  mvaMETToken_(consumes<pat::METCollection>(mvaMETTag_))
  mvaMETToken_(consumes<reco::PFMETCollection>(mvaMETTag_))
{
  produces<std::vector<vhtm::MET>>("vhtmPFMETVector").setBranchAlias("vhtmPFMETVector");
  produces<std::vector<vhtm::MET>>("vhtmCorrMETVector").setBranchAlias("vhtmCorrMETVector");
  produces<std::vector<vhtm::MET>>("vhtmPuppiMETVector").setBranchAlias("vhtmPuppiMETVector");
  produces<std::vector<vhtm::MET>>("vhtmMVAMETVector").setBranchAlias("vhtmMVAMETVector");
}
void METBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  pfList_ = new std::vector<vhtm::MET>();
  tree->Branch("MET", "std::vector<vhtm::MET>", &pfList_, 32000, -1);
  tree->Branch("nMET", &fnPFMET_, "fnPFMET_/I");

  corrList_ = new std::vector<vhtm::MET>();
  tree->Branch("corrMET", "std::vector<vhtm::MET>", &corrList_, 32000, -1);
  tree->Branch("corrnMET", &fnCorrMET_, "fnCorrMET_/I");

  puppiList_ = new std::vector<vhtm::MET>();
  tree->Branch("puppiMET", "std::vector<vhtm::MET>", &puppiList_, 32000, -1);
  tree->Branch("puppinMET", &fnPuppiMET_, "fnPuppiMET_/I");

  mvaList_ = new std::vector<vhtm::MET>();
  tree->Branch("mvaMET", "std::vector<vhtm::MET>", &mvaList_, 32000, -1);
  tree->Branch("mvanMET", &fnMVAMET_, "fnMVAMET_/I");
}
void METBlock::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  fillMET(iEvent, iSetup, pfList_, fnPFMET_, pfMETTag_, pfMETToken_);
  fillMET(iEvent, iSetup, corrList_, fnCorrMET_, corrMETTag_, corrMETToken_);
  fillMET(iEvent, iSetup, puppiList_, fnPuppiMET_, puppiMETTag_, puppiMETToken_);
  fillMET(iEvent, iSetup, mvaList_, fnMVAMET_, mvaMETTag_, mvaMETToken_);
  //put the vhtm collections in edm
  std::auto_ptr<std::vector<vhtm::MET>> pv1(new std::vector<vhtm::MET>(*pfList_));
  iEvent.put(pv1,"vhtmPFMETVector");
  std::auto_ptr<std::vector<vhtm::MET>> pv2(new std::vector<vhtm::MET>(*corrList_));
  iEvent.put(pv2,"vhtmCorrMETVector");
  std::auto_ptr<std::vector<vhtm::MET>> pv3(new std::vector<vhtm::MET>(*puppiList_));
  iEvent.put(pv3,"vhtmPuppiMETVector");
  std::auto_ptr<std::vector<vhtm::MET>> pv4(new std::vector<vhtm::MET>(*mvaList_));
  iEvent.put(pv4,"vhtmMVAMETVector");
}
void METBlock::fillMET(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       std::vector<vhtm::MET>* list,
                       int& nMET,
                       const edm::InputTag& iTag,
                       const edm::EDGetTokenT<pat::METCollection>& token)
{
  // Reset the TClonesArray and the nObj variables
  list->clear();
  nMET = 0;

  edm::Handle<pat::METCollection> metColl;
  bool found = iEvent.getByToken(token, metColl);

  if (found && metColl.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << metColl->size();
    for (const pat::MET& v: *metColl) {
      if (list->size() == kMaxMET_) {
        edm::LogInfo("METBlock") << "Too many PAT MET, nMET = " << list->size()
				 << ", label: " << iTag;
        break;
      }
      vhtm::MET mobj;
      // fill in all the vectors
      mobj.met          = v.pt();
      mobj.metphi       = v.phi();
      mobj.sumet        = v.sumEt();
      mobj.metuncorr    = v.uncorPt();
      mobj.metphiuncorr = v.uncorPhi();
      //mobj.sumetuncorr  = v.sumEt() - v.corSumEt(pat::MET::uncorrALL);

      list->push_back(mobj);
    }
    nMET = list->size();      
  }
  else {
    edm::LogError("METBlock") << "Error! Failed to get pat::MET collection for label: "
                              << iTag;
  }
}
void METBlock::fillMET(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       std::vector<vhtm::MET>* list,
                       int& nMET,
                       const edm::InputTag& iTag,
                       const edm::EDGetTokenT<reco::PFMETCollection>& token)
{
  // Reset the TClonesArray and the nObj variables
  list->clear();
  nMET = 0;

  edm::Handle<reco::PFMETCollection> metColl;
  bool found = iEvent.getByToken(token, metColl);

  if (found && metColl.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << metColl->size();
    for (const reco::PFMET& v: *metColl) {
      if (list->size() == kMaxMET_) {
        edm::LogInfo("METBlock") << "Too many PAT MET, nMET = " << list->size()
				 << ", label: " << iTag;
        break;
      }
      vhtm::MET mobj;
      // fill in all the vectors
      mobj.met          = v.pt();
      mobj.metphi       = v.phi();
      mobj.sumet        = v.sumEt();
      //mobj.metuncorr    = v.uncorPt();
      //mobj.metphiuncorr = v.uncorPhi();
      //mobj.sumetuncorr  = v.sumEt() - v.corSumEt(reco::PFMET::uncorrALL);

      list->push_back(mobj);
    }
    nMET = list->size();      
  }
  else {
    edm::LogError("METBlock") << "Error! Failed to get pat::MET collection for label: "
                              << iTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);

