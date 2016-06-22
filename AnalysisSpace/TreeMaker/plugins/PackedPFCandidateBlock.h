#ifndef __AnalysisSpace_TreeMaker_PackedPFCandidateBlock_h
#define __AnalysisSpace_TreeMaker_PackedPFCandidateBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

class PackedPFCandidateBlock: public edm::one::EDAnalyzer<edm::one::SharedResources> 
{
private:
  virtual void beginJob();
  //virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

public:
  explicit PackedPFCandidateBlock(const edm::ParameterSet& iConfig);
  virtual ~PackedPFCandidateBlock();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void calcIsoFromPF(const pat::PackedCandidate& v,
                                           edm::Handle<pat::PackedCandidateCollection>& pfs,
                                           double cone, std::vector<double>& iso);
  enum {
    kMaxPackedPFCandidate = 250
  };

private:
  const int verbosity_;
  std::vector<vhtm::PackedPFCandidate>* list_;
  int fnPackedPFCandidate_;

  const edm::InputTag pfcandTag_;
  std::vector<int> pdgTosave_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

};
#endif
