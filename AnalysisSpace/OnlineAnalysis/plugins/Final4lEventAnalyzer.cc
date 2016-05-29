#include <memory>
#include<iomanip>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "AnalysisSpace/OnlineAnalysis/plugins/Final4lEventAnalyzer.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

#include "ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h"

//header file for tree
#include "AnalysisSpace/TreeMaker/interface/Utility.h"



using namespace MEMNames;
//
// constructors and destructor
//
Final4lEventAnalyzer::Final4lEventAnalyzer(const edm::ParameterSet& iConfig) :
  ZZcandTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZZcandSrc")),
  lepCleanedlooseJetTag_(iConfig.getUntrackedParameter<edm::InputTag>("lepCleanedJetSrc")),
  tightIsoEleFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightIsoElectronFSRSrc")),
  tightIsoMuFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightIsoMuonFSRSrc")),
  ZZcandToken_(consumes<std::vector<vhtm::ZZcandidate>>(ZZcandTag_)),
  lepCleanedlooseJetToken_(consumes<pat::JetCollection>(lepCleanedlooseJetTag_)),
  tightIsoElectronFSRToken_(consumes<std::vector<vhtm::IsoElectron>>(tightIsoEleFSRpairTag_)),
  tightIsoMuonFSRToken_(consumes<std::vector<vhtm::IsoMuon>>(tightIsoMuFSRpairTag_)),
  syncFilename_(iConfig.getUntrackedParameter<std::string>("syncFilename")),
  runStandalone_(iConfig.getUntrackedParameter<bool>("runStandalone")),
  computeKD_(iConfig.getUntrackedParameter<bool>("computeKD",false)),
  isData_(iConfig.getUntrackedParameter<bool>("isData",false))
{
   syncDumpf_.open(syncFilename_.c_str(), ios::out);
   if (!syncDumpf_) {
    std::cerr << "Output File: " << syncFilename_ << " could not be opened!" << endl;
   }
   std::cout << "Is Final4l Analyzer running in standalone mode?>>" << runStandalone_ 
             << "\tcomputeKD?>>" << computeKD_
             << "\tisData?>>" << isData_ 
             << std::endl;
   //now do what ever initialization is needed
   //usesResource("TFileService");
}


Final4lEventAnalyzer::~Final4lEventAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete selev_;
   
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Final4lEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //std::cout << "DEBUG P1" << std::endl;
   selev_->reset();
   //using namespace edm;
   edm::Handle<pat::JetCollection> jetVec;
   //bool jetFound = 
   iEvent.getByToken(lepCleanedlooseJetToken_, jetVec);
   int nJets = jetVec->size();
   //std::cout << "Jet size=" << nJets << std::endl;
   
   edm::Handle<std::vector<vhtm::ZZcandidate>> ZZcandColl;
   //bool zzFound = 
   iEvent.getByToken(ZZcandToken_, ZZcandColl);
   //std::cout << "ZZcandColl size=" << ZZcandColl->size() << std::endl;

   edm::Handle<std::vector<vhtm::IsoElectron>> tightIsoEle;
   iEvent.getByToken(tightIsoElectronFSRToken_, tightIsoEle);

   edm::Handle<std::vector<vhtm::IsoMuon>> tightIsoMu;
   iEvent.getByToken(tightIsoMuonFSRToken_, tightIsoMu);

   //std::cout << "DEBUG P2" << std::endl;
   
   if(!ZZcandColl->empty()) {
     vhtm::ZZcandidate selectedZZ;
     selectBestZZCandidate(ZZcandColl,selectedZZ);

   std::cout << "DEBUG P3" << std::endl;
     TLorentzVector jet1P4, jet2P4;
     jet1P4.SetPtEtaPhiE(0.,0.,0.,0.);
     jet2P4.SetPtEtaPhiE(0.,0.,0.,0.);
     if (nJets) {
       jet1P4 = HZZ4lUtil::getP4(jetVec->at(0));
       if (nJets > 1) {
         jet2P4 = HZZ4lUtil::getP4(jetVec->at(1));
       }
     }
     int nbJets = 0;
     for(auto& j : *jetVec)
       if(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.800)   nbJets++;
     if(computeKD_)   computeKD(selectedZZ, jet1P4, jet2P4, nJets, selev_->kd);

   //std::cout << "DEBUG P4" << std::endl;
     double m4lrefit = -1.;
     double m4lrefiterr = -1.;
     doZmassrefit(selectedZZ, m4lrefit, m4lrefiterr);
     int nextlep = findExtraLeptons(*tightIsoEle, *tightIsoMu, selectedZZ);
     std::cout << "Extra Lepton=" << nextlep << std::endl;
     //compute this only if nlepton = 4 and if #jets >= 2  
     bool jetPaircond = false;
     if(nextlep == 0 && jetVec->size() >= 2) jetPaircond = hasJetPair(*jetVec);
   //std::cout << "DEBUG P5" << std::endl;
     //compute EventCategory
     TLorentzVector final4lP4;
     for(auto& p : selectedZZ.lepP4Vec )  final4lP4+=p;
     for(auto& p : selectedZZ.fsrP4Vec )  final4lP4+=p;
     int cat = findEventCategory(nextlep+4, jet1P4, jet2P4, final4lP4, jetVec->size(), nbJets, jetPaircond);
   //std::cout << "DEBUG P6" << std::endl;
     //dump file 
     HZZ4lUtil::syncDumper(iEvent.id().run(), iEvent.id().luminosityBlock(),iEvent.id().event(), selectedZZ, 
                           nJets, jet1P4.Pt(), jet2P4.Pt(), selev_->kd, cat, m4lrefit, m4lrefiterr, 0,syncDumpf_);
     //if(!runStandalone_) {
       selev_->run = iEvent.id().run();
       selev_->lumi = iEvent.id().luminosityBlock();
       selev_->event = iEvent.id().event();
       selev_->mass4l = selectedZZ.m4l;
       selev_->mZ1 = selectedZZ.mZ1;
       selev_->mZ2 = selectedZZ.mZ2;
       //selev_->kd;//D_bkg^kin,D_bkg,D_gg,D_HJJ^VBF,D_0
       selev_->nJets = nJets;
       selev_->jet1pt = jet1P4.Pt();
       selev_->jet2pt = jet2P4.Pt();
       //selev_->category;
       //selev_->m4lRefit;
       //selev_->m4lRefitError;
       //selev_->weight;
     //}
   }
   //std::cout << "DEBUG P6" << std::endl;
}
//If more than one ZZ candidate survives, choose the one with the Z1 closest in mass to nominal Z. 
//If two or more candidates include the same Z1 and differ just by the Z2, choose the one with the highest-pT Z2 leptons. 
//(Other Z pairing choices, e.g. by best Dbkg, could be considered later). 
void 
Final4lEventAnalyzer::selectBestZZCandidate(const edm::Handle<std::vector<vhtm::ZZcandidate> >& zzlist, vhtm::ZZcandidate& Zcand) {
  Zcand = zzlist->at(0);
  for(unsigned int i = 1; i<zzlist->size(); i++) {
    if(zzlist->at(i).Z1mdiff < Zcand.Z1mdiff) 
      Zcand = zzlist->at(i);
    else if(zzlist->at(i).Z1mdiff == Zcand.Z1mdiff) {
      if(zzlist->at(i).Z2lepPtavg > Zcand.Z2lepPtavg)    Zcand = zzlist->at(i);
    }
  } 
}
//function to compute the Kinematic Discriminants
void Final4lEventAnalyzer::computeKD(const vhtm::ZZcandidate& ZZcand, const TLorentzVector& jet1P4, const TLorentzVector& jet2P4,
			             int nJets,std::map<std::string,double>&  kd) 
{
  std::vector<TLorentzVector> partP = ZZcand.lepP4Vec;
  std::vector<int> partId = ZZcand.lepcharge;
  std::vector<TLorentzVector> partPprod = ZZcand.lepP4Vec;
  partPprod.push_back(jet1P4);
  partPprod.push_back(jet2P4); // Can also use partPprod.push_back(nullFourVector) instead for integrated VBF MEs
  //TLorentzVector nullFourVector(0, 0, 0, 0);
  
  std::vector<int> partIdprod = ZZcand.lepcharge; 
  partIdprod.push_back(0); // For leptonic ZH in the future, this could actually be another lepton flavor
  partIdprod.push_back(0); // For leptonic ZH in the future, this could actually be the opposite lepton flavor
  
  double p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
  double p0minus_VAJHU, Dgg10_VAMCFM;
  double phjj_VAJHU, pvbf_VAJHU;
  
  MEMs combinedMEM(13, 125, "CTEQ6L");
  combinedMEM.computeME(MEMNames::kSMHiggs, MEMNames::kJHUGen, partP, partId, p0plus_VAJHU); // Calculation of SM gg->H->4l JHUGen ME
  combinedMEM.computeME(MEMNames::k0minus, MEMNames::kJHUGen, partP, partId, p0minus_VAJHU); // Calculation of PS (0-, fa3=1) gg->H->4l JHUGen ME 
  combinedMEM.computeME(MEMNames::kggHZZ_10, MEMNames::kMCFM, partP, partId, Dgg10_VAMCFM); // Direct calc of Dgg (D^kin for off-shell) from MCFM MEs
  combinedMEM.computeME(MEMNames::kqqZZ, MEMNames::kMCFM, partP, partId, bkg_VAMCFM); // qq->4l background calculation from MCFM
  combinedMEM.computePm4l(partP,partId, MEMNames::kNone, p0plus_m4l, bkg_m4l); // m4l probabilities for signal and background, nominal resolution
  combinedMEM.computeME(MEMNames::kJJ_SMHiggs_GG, MEMNames::kJHUGen, partPprod, partIdprod, phjj_VAJHU); // SM gg->H+2j
  combinedMEM.computeME(MEMNames::kJJ_SMHiggs_VBF, MEMNames::kJHUGen, partPprod, partIdprod, pvbf_VAJHU);  // SM VBF->H
  
  // Dgg already obtained above
  float D_bkg_kin = p0plus_VAJHU / ( p0plus_VAJHU + bkg_VAMCFM ); // D^kin_bkg
  float D_bkg = p0plus_VAJHU * p0plus_m4l / ( p0plus_VAJHU * p0plus_m4l + bkg_VAMCFM * bkg_m4l ); // D^kin including superMELA
  float D_g4 = p0plus_VAJHU / ( p0plus_VAJHU + p0minus_VAJHU ); // D_0-
  float Djet_VAJHU = pvbf_VAJHU / ( pvbf_VAJHU + phjj_VAJHU ); // D^VBF_HJJ
  
  kd["Dgg10_VAMCFM"] = Dgg10_VAMCFM; 
  kd["D_bkg_kin"] = D_bkg_kin; 
  kd["D_bkg"] = D_bkg; 
  kd["D_g4"] = D_g4; 
  kd["Djet_VAJHU"] = (nJets >= 2) ? Djet_VAJHU : -1.;
}
//Function to do mass refit
//implementation https://github.com/tocheng/KinZfitter
void Final4lEventAnalyzer::doZmassrefit(vhtm::ZZcandidate& ZZcand, double& mass4lREFIT , double& mass4lErrREFIT) 
{
  vector<reco::Candidate *> selectedLeptons;
  TLorentzVector pL11, pL12, pL21, pL22;
  std::map<unsigned int, TLorentzVector> selectedFsrMap;
  TLorentzVector tempP4;
  tempP4.SetPtEtaPhiE(0.,0.,0.,0.);
  if(ZZcand.flavour == HZZ4lUtil::ZZType::eeee) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep2));
  } else if(ZZcand.flavour == HZZ4lUtil::ZZType::mmmm) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep2));
    if(ZZcand.Z2mm.lep2hasfsr)   selectedFsrMap[3] = HZZ4lUtil::getP4(ZZcand.Z2mm.fsrl2);   
  } else if(ZZcand.flavour == HZZ4lUtil::ZZType::mmee) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep2));
  } else if(ZZcand.flavour == HZZ4lUtil::ZZType::eemm) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep2));
  }
  selectedFsrMap[0] = ZZcand.fsrP4Vec[0];
  selectedFsrMap[1] = ZZcand.fsrP4Vec[1];
  selectedFsrMap[2] = ZZcand.fsrP4Vec[2];
  selectedFsrMap[3] = ZZcand.fsrP4Vec[3];
  kinZfitter_->Setup(selectedLeptons, selectedFsrMap);
  kinZfitter_->KinRefitZ1();
  // refit mass4l
  mass4lREFIT = kinZfitter_->GetRefitM4l();
  // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
  //vector < TLorentzVector > p4 = kinZfitter->GetRefitP4s();//use later 
  // refitted mass4l error
  mass4lErrREFIT = kinZfitter_->GetRefitM4lErrFullCov();
}
// Function to calculate number of extra tight leptons passing isolation apart from leptons in
// ZZ candidate
int Final4lEventAnalyzer::findExtraLeptons(const std::vector<vhtm::IsoElectron>& ieVec, const std::vector<vhtm::IsoMuon>& imVec, const vhtm::ZZcandidate& ZZcand)
{
  int nExtralep = 0;
  const TLorentzVector& pL11 = ZZcand.lepP4Vec[0]; 
  const TLorentzVector& pL12 = ZZcand.lepP4Vec[1]; 
  const TLorentzVector& pL21 = ZZcand.lepP4Vec[2]; 
  const TLorentzVector& pL22 = ZZcand.lepP4Vec[3];

  for(auto& ie:ieVec) {
    TLorentzVector eleP4 = HZZ4lUtil::getP4(ie.ele);
    if (eleP4 == pL11 || eleP4 == pL12 || eleP4 == pL21 || eleP4 == pL22) continue;
    nExtralep++;
  }
  for(auto& im:imVec) {
    TLorentzVector muP4 = HZZ4lUtil::getP4(im.mu);
    if (muP4 == pL11 || muP4 == pL12 || muP4 == pL21 || muP4 == pL22) continue;
    nExtralep++;
  }
  return nExtralep;
}
//Find event category
/*
 exactly 4 leptons + at least 2 jets + at most 1 b-tag jet in the event +Djet>0.5
    ➔ VBF-tagged category (#2)
 exactly 4 leptons + at least one pair of jets passing { |η|<2.4, pT>40, 60<mjj<120 } + pT,4l>m4l,
    OR: exactly 4 leptons + exactly 2 jets, both with b-tag
    ➔ VH-hadronic-tagged category (#4)
 no more than 2 jets + no b-tagged jet in the event; + at least 5 leptons
    ➔ VH-leptonic-tagged category (#3)
 at least 3 jets, of which at least 1 b-tagged,
    OR: at least 5 leptons
    ➔ ttH-tagged category (#5)
 at least 1 jet
    ➔ 1-jet-tagged category (#1)
 remaining events
    ➔ untagged category (#0) 
*/
int Final4lEventAnalyzer::findEventCategory(const int nleptons, const TLorentzVector& j1P4, const TLorentzVector& j2P4, 
                                            const TLorentzVector& m4lP4withfsr, const int njets, const int nbjets, const bool jetPaircond)
{
  double djet = -1., mjj = 0.;
  if (njets >= 2) {
    mjj = (j1P4 + j2P4).M();
    djet = 0.18 * std::fabs(j1P4.Eta() - j2P4.Eta()) + 1.92E-04 * mjj;
  }
  int cat = 0;
  if (nleptons == 4 && njets >= 2 && nbjets <= 1 && djet > 0.5) 
    cat = 2;
  else if ( (nleptons == 4 && 
	     njets >= 2 && jetPaircond &&
	     m4lP4withfsr.Pt() > m4lP4withfsr.M() ) ||
            (nleptons == 4 && njets == 2 && nbjets == 2) )
    cat = 4;
  else if (njets <= 2 && nbjets == 0 && nleptons >= 5)
    cat = 3;
  else if ( (njets >= 3 && nbjets >= 1) || nleptons >= 5 )
    cat = 5;
  else if (njets >= 1)
    cat = 1;
  return cat;
}

bool Final4lEventAnalyzer::hasJetPair(const pat::JetCollection& jetList) {
  for (unsigned int i = 0; i < jetList.size(); ++i) {
    auto const& j1 = jetList[i];
    auto const& j1P4 = HZZ4lUtil::getP4(j1);
    for (unsigned int j = i+1; j < jetList.size(); ++j) {
      auto const& j2 = jetList[j];
      auto const& j2P4 = HZZ4lUtil::getP4(j2);
      double mjj = (j1P4+j2P4).M();      
      if ( (std::fabs(j1P4.Eta()) < 2.4 && j1P4.Pt() > 40.) && 
	   (std::fabs(j2P4.Eta()) < 2.4 && j2P4.Pt() > 40.) && 
	   (mjj > 60. && mjj < 120.) ) return true;
    }
  }
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Final4lEventAnalyzer::beginJob()
{
  selev_ = new vhtm::SelectedEvent();
  if(!runStandalone_) {
    TTree* tree = vhtm::Utility::getTree("vhtree");
    tree->Branch("SelectedEvent", "vhtm::SelectedEvent", &selev_, 32000, -1);
    std::cout << "Tree pointer found>>" << tree << "\t branch booked" << std::endl;
  }
  kinZfitter_ = new KinZfitter(isData_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Final4lEventAnalyzer::endJob() 
{
  delete kinZfitter_;
  syncDumpf_.close();
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Final4lEventAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Final4lEventAnalyzer);
