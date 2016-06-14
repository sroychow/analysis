#ifndef __AnalysisSpace_OnlineAnalysis_ZTnpAnalyzer_h
#define __AnalysisSpace_OnlineAnalysis_ZTnpAnalyzer_h
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


#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
using EleFSRCollection=std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>;
using MuFSRCollection=std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>;

class ZTnpAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZTnpAnalyzer(const edm::ParameterSet&);
      ~ZTnpAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void fillZeeTnp(const edm::Handle<EleFSRCollection>& looseSIPEleFSRpair, const edm::Handle<pat::PackedCandidateCollection>& fsr,const reco::Vertex& vit, 
                      const double& evrho);
      void fillZmumuTnp(const edm::Handle<MuFSRCollection>& looseSIPMuFSRpair, const edm::Handle<pat::PackedCandidateCollection>& fsr,const reco::Vertex& vit);

      //virtual void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& iSetup) override;
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      int verbosity_;
      const edm::InputTag fsrTag_;
      const edm::InputTag looseSIPEleFSRpairTag_;
      const edm::InputTag looseSIPMuFSRpairTag_;
      const edm::InputTag fixedGridRhoFastjetAllTag_;
      const edm::InputTag vertexTag_;

      const edm::EDGetTokenT<pat::PackedCandidateCollection> fsrToken_;
      const edm::EDGetTokenT<EleFSRCollection> looseElectronFSRToken_;
      const edm::EDGetTokenT<MuFSRCollection> looseMuonFSRToken_;
      const edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
      const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      std::vector<vhtm::ZtnP>*   zlltnp_;
      bool isData_;
};
#endif
