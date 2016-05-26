#ifndef __AnalysisSpace_OnlineAnalysis_FSRBlock_h
#define __AnalysisSpace_OnlineAnalysis_FSRBlock_h
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
//Root
#include "TLorentzVector.h"

//
// class declaration
//

class FSRBlock : public edm::stream::EDProducer<> {
   public:
      explicit FSRBlock(const edm::ParameterSet&);
      ~FSRBlock();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      
      void calcIsoFromPF(const pat::PackedCandidate& v, 
                         edm::Handle<pat::PackedCandidateCollection>& pfs, 
                         double cone, std::vector<double>& iso);
      bool passedSuperClusterVeto(const pat::PackedCandidate& pfcand, edm::Handle<pat::ElectronCollection>& le, bool verbose=false);
      double findClosestLepton(const pat::PackedCandidate& pfPho, edm::Handle<pat::ElectronCollection>& le, 
                               edm::Handle<pat::MuonCollection>& lmu,int& muindx, int& elindx);
      template<class T>
      TLorentzVector getP4(const T& pf) {
        TLorentzVector lv;
        lv.SetPtEtaPhiE(pf.pt(), pf.eta(), pf.phi(), pf.energy());
        return lv;
      }

      int fnPhoton_;

      int verbosity_;
      const edm::InputTag pfcandTag_;
      const edm::InputTag looseElectronTag_;
      const edm::InputTag looseMuonTag_;
      const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      const edm::EDGetTokenT<pat::ElectronCollection> looseElectronToken_;
      const edm::EDGetTokenT<pat::MuonCollection> looseMuonToken_;
      const std::string selectedFSRColl_;
      const std::string looseSIPEleFSRColl_;
      const std::string looseSIPMuFSRColl_;
      std::vector<pat::PackedCandidate>* FSRVec_;
      std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>* looseEleFSRpairVec_;
      std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>*     looseMuFSRpairVec_;
};
#endif
