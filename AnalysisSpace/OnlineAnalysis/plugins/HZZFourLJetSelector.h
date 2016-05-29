// system include files
#ifndef __AnalysisSpace_OnlineAnalysis_HZZFourLJetSelector_h
#define __AnalysisSpace_OnlineAnalysis_HZZFourLJetSelector_h
#include <memory>
#include <string>
#include <vector>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"

//
// class declaration
//

class HZZFourLJetSelector : public edm::stream::EDProducer<> {
   public:
      explicit HZZFourLJetSelector(const edm::ParameterSet&);
      ~HZZFourLJetSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob(/*edm::StreamID*/) /*override*/;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() /*override*/;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
      void calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Muon& v, std::vector<double>& iso);     
      bool jetLeptonCleaning(const TLorentzVector& jet, const std::vector<TLorentzVector>& lepfsrp4vec,double dR);
      bool isLooseJet(const pat::Jet& jet); 
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

   private:

      //const int verbosity_;

      const edm::InputTag tightIsoEleFSRpairTag_;
      const edm::InputTag tightIsoMuFSRpairTag_;
      const edm::InputTag jetTag_;

      const edm::EDGetTokenT<std::vector<vhtm::IsoElectron>> tightIsoElectronFSRToken_;
      const edm::EDGetTokenT<std::vector<vhtm::IsoMuon>> tightIsoMuonFSRToken_;
      const edm::EDGetTokenT<pat::JetCollection> jetToken_;
   
      const std::string looseJetOutputColl_;
      std::vector<pat::Jet>* looseJetVec_;
};
#endif
