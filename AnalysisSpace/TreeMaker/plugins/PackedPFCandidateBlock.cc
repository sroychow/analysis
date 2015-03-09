#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "AnalysisSpace/TreeMaker/plugins/PackedPFCandidateBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

// Constructor
PackedPFCandidateBlock::PackedPFCandidateBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfcandTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfCands",edm::InputTag("packedPFCandidates"))), 
  pdgTosave_(iConfig.getParameter<std::vector<int>>("pdgTosave")),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_))
{}
void PackedPFCandidateBlock::beginJob() 
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::PackedPFCandidate>();
  tree->Branch("PackedPFCandidate", "std::vector<vhtm::PackedPFCandidate>", &list_, 32000, 2);
  tree->Branch("nPackedPFCandidate", &fnPackedPFCandidate_, "fnPackedPFCandidate_/I");
}
void PackedPFCandidateBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  list_->clear();
  fnPackedPFCandidate_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);

  
  if ( pfs.isValid()) {
    edm::LogInfo("PackedPFCandidateBlock") << "Total # PackedPFCandidate: " << pfs->size();
    for (pat::PackedCandidate const& v: *pfs) {
      //if (list_->size() == kMaxPhoton) {
      //	edm::LogInfo("PackedPFCandidate") << "Too many PAT Photon, fnPhoton = " 
      //                              << fnPhoton_; 
      //	break;
      //}

      int pdg = std::abs(v.pdgId());
      if( std::find(pdgTosave_.begin(),  pdgTosave_.end(), pdg) == pdgTosave_.end() ) continue;
      
      //if( pdg != std::abs(11) && pdg != std::abs(13) &&
      //    pdg != std::abs(15) && pdg != std::abs(22) 
      //  ) continue;

      vhtm::PackedPFCandidate pfCand;
      
      pfCand.pt = v.pt();
      pfCand.eta = v.eta();
      pfCand.phi = v.phi();
      pfCand.energy = v.energy();
      
      pfCand.pdgId = v.pdgId();
      pfCand.charge = v.charge();
   
      pfCand.vx = v.vx();
      pfCand.vz = v.vz();
      pfCand.vz = v.vz();
     
      pfCand.fromPV = v.fromPV();
      pfCand.dxy = v.dxy();
      pfCand.dz = v.dz();
      pfCand.dxyError = v.dxyError();   
      pfCand.dzError = v.dzError();   
   
      list_->push_back(pfCand);
    }
    fnPackedPFCandidate_ = list_->size();
  }
  else {
    edm::LogError("PackedPFCandidateBlock") << "Error >> Failed to get pat::Photon for label: " 
                                 << pfcandTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PackedPFCandidateBlock);

