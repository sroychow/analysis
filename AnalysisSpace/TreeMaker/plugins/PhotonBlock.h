#ifndef __AnalysisSpace_TreeMaker_PhotonBlock_h
#define __AnalysisSpace_TreeMaker_PhotonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

namespace vhtm {
  class Photon;
}

class PhotonBlock : public edm::EDAnalyzer 
{
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

 public:
  explicit PhotonBlock(const edm::ParameterSet& iConfig);
  virtual ~PhotonBlock() {}

  enum {
    kMaxPhoton = 100
  };
 private:
  std::vector<vhtm::Photon>* list_;
  int fnPhoton_;

  int verbosity_;
  const edm::InputTag photonTag_;
  const edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
};
#endif
