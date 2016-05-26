#ifndef __AnalysisSpace_OnlineAnalysis_ZCandidateProducer_h
#define __AnalysisSpace_OnlineAnalysis_ZCandidateProducer_h
// system include files
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
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"


#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"

//Root
#include "TLorentzVector.h"

//
// class declaration
//

class ZCandidateProducer : public edm::stream::EDProducer<> {
   public:
      explicit ZCandidateProducer(const edm::ParameterSet&);
      ~ZCandidateProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      void ZmumuSelector(const std::vector<vhtm::IsoMuon>& lepPhotonPairVec, std::vector<vhtm::Zmumu>& candList);
      void ZeeSelector(const std::vector<vhtm::IsoElectron>& lepPhotonPairVec, std::vector<vhtm::Zee>& candList);


      int verbosity_;
      const edm::InputTag tightIsoEleFSRpairTag_;
      const edm::InputTag tightIsoMuFSRpairTag_;
      

      const edm::EDGetTokenT<std::vector<vhtm::IsoElectron>> tightIsoElectronFSRToken_;
      const edm::EDGetTokenT<std::vector<vhtm::IsoMuon>> tightIsoMuonFSRToken_;

      const std::string ZeeColl_;
      const std::string ZmumuColl_;
     
      std::vector<vhtm::Zee>* ZeeVec_;
      std::vector<vhtm::Zmumu>* ZmumuVec_;
};
#endif
