// system include files
#ifndef __AnalysisSpace_OnlineAnalysis_HZZFourLMuonSelector_h
#define __AnalysisSpace_OnlineAnalysis_HZZFourLMuonSelector_h
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"


//
// class declaration
//

class HZZFourLMuonSelector : public edm::stream::EDProducer<> {
   public:
      explicit HZZFourLMuonSelector(const edm::ParameterSet&);
      ~HZZFourLMuonSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob(/*edm::StreamID*/) /*override*/;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() /*override*/;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
      bool isTrackerHighPt(const pat::Muon & mu, const reco::Vertex & primaryVertex);
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
   public:

   private:
      std::vector<pat::Muon>* looselist_;
      std::vector<pat::Muon>* looseSIPlist_;
      std::vector<pat::Muon>* tightlist_;
      int fnMuon_;

      const int verbosity_;

      const edm::InputTag muonTag_;
      const edm::InputTag vertexTag_;

      const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

      const std::string looseMuOutputColl_;
      const std::string looseSIPMuOutputColl_;
      const std::string tightMuOutputColl_;
};
#endif
