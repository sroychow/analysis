#ifndef __AnalysisSpace_OnlineAnalysis_Final4lEventAnalyzer_h
#define __AnalysisSpace_OnlineAnalysis_Final4lEventAnalyzer_h
// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"

//header for kinematic refit
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Final4lEventAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Final4lEventAnalyzer(const edm::ParameterSet&);
      ~Final4lEventAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void selectBestZZCandidate(const edm::Handle<std::vector<vhtm::ZZcandidate> >& zzlist, vhtm::ZZcandidate& Zcand);
      void computeKD(const vhtm::ZZcandidate& ZZcand, const TLorentzVector& jet1P4, const TLorentzVector& jet2P4,
	             int nJets,std::map<std::string,double>&  kd);
      void doZmassrefit(vhtm::ZZcandidate& ZZcand, double& mass4lREFIT, double& mass4lErrREFIT);
      int findExtraLeptons(const std::vector<vhtm::IsoElectron>& ieVec, const std::vector<vhtm::IsoMuon>& imVec, const vhtm::ZZcandidate& ZZcand);
      int findEventCategory(const int nleptons, const TLorentzVector& j1P4, const TLorentzVector& j2P4, 
                            const TLorentzVector& m4lP4withfsr, const int njets, const int nbjets, const bool jetPaircond);
      bool hasJetPair(const pat::JetCollection& jetList);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) override;
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      const edm::InputTag ZZcandTag_;
      const edm::InputTag lepCleanedlooseJetTag_;
      const edm::InputTag tightIsoEleFSRpairTag_;
      const edm::InputTag tightIsoMuFSRpairTag_;

      const edm::EDGetTokenT<std::vector<vhtm::ZZcandidate>> ZZcandToken_;
      const edm::EDGetTokenT<pat::JetCollection> lepCleanedlooseJetToken_;
      const edm::EDGetTokenT<std::vector<vhtm::IsoElectron>> tightIsoElectronFSRToken_;
      const edm::EDGetTokenT<std::vector<vhtm::IsoMuon>> tightIsoMuonFSRToken_;

      const std::string syncFilename_;
      std::ofstream syncDumpf_;
      vhtm::SelectedEvent* selev_;
      KinZfitter *kinZfitter_;
      bool runStandalone_; 
      bool computeKD_;
      bool isData_;
};
#endif
