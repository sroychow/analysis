#ifndef __AnalysisSpace_OnlineAnalysis_FourLeptonSelector_h
#define __AnalysisSpace_OnlineAnalysis_FourLeptonSelector_h

/*
\class FourLeptonSelector FourLeptonSelector.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Suvankar Roy Chowdhury
//         Created:  Mon, 06 Feb 2017 11:56:27 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//header for kinematic refit
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"

//header for mela
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "AnalysisSpace/OnlineAnalysis/src/KDCalculator.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class FourLeptonSelector : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit FourLeptonSelector(const edm::ParameterSet&);
      ~FourLeptonSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      bool verbosity_;
      bool isMC_;
      const std::string syncFilename_;
      std::ofstream syncDumpf_;
      int nEvtswith4Lep_;
      int nEvtwith2Z_;
      int nEvtwithZZ_;
      int nEvtwith2e2m_;
      int nEvtwith4m_;
      int nEvtwith4e_;
      const edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
      const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      //Selected Electrons token
      const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      const edm::EDGetTokenT<pat::JetCollection> jetToken_;
      Mela* mela_;
      //variables for saving in tree if ZZ candidate is found
      unsigned long int run_;
      unsigned long int lumi_;
      unsigned long int event_;
      bool hasSelectedZZ_;
      float m4l_;
      float mZ1_;
      float mZ2_;
};
#endif
