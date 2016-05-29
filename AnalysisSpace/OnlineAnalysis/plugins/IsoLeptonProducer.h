#ifndef __AnalysisSpace_OnlineAnalysis_IsoLeptonProducer_h
#define __AnalysisSpace_OnlineAnalysis_IsoLeptonProducer_h
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

#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//Root
#include "TLorentzVector.h"

//
// class declaration
//

using EleFSRCollection=std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>;
using MuFSRCollection=std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>;
class IsoLeptonProducer : public edm::stream::EDProducer<> {
   public:
      explicit IsoLeptonProducer(const edm::ParameterSet&);
      ~IsoLeptonProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      
      template<class T>
      TLorentzVector getP4(const T& pf) {
        TLorentzVector lv;
        lv.SetPtEtaPhiE(pf.pt(), pf.eta(), pf.phi(), pf.energy());
        return lv;
      }

      int verbosity_;
      const edm::InputTag fsrTag_;
      const edm::InputTag looseSIPEleFSRpairTag_;
      const edm::InputTag looseSIPMuFSRpairTag_;
      const edm::InputTag fixedGridRhoFastjetAllTag_;

      const edm::EDGetTokenT<pat::PackedCandidateCollection> fsrToken_;
      const edm::EDGetTokenT<EleFSRCollection> looseElectronFSRToken_;
      const edm::EDGetTokenT<MuFSRCollection> looseMuonFSRToken_;
      const edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;

      const std::string tightSIPEleFSRColl_;
      const std::string tightSIPMuFSRColl_;
      std::vector<vhtm::IsoElectron>* tightEleFSRpairVec_;
      std::vector<vhtm::IsoMuon>* tightMuFSRpairVec_;
};
#endif
