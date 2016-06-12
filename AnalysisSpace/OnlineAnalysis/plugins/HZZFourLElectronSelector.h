#ifndef __AnalysisSpace_OnlineAnalysis_HZZFourLElectronSelector_h
#define __AnalysisSpace_OnlineAnalysis_HZZFourLElectronSelector_h
// -*- C++ -*-
//
// Package:    HZZFourLElectronSelector/HZZFourLElectronSelector
// Class:      HZZFourLElectronSelector
// 
/**\class HZZFourLElectronSelector HZZFourLElectronSelector.cc HZZFourLElectronSelector/HZZFourLElectronSelector/plugins/HZZFourLElectronSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


//
// class declaration
//
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

class HZZFourLElectronSelector : public edm::stream::EDProducer<> {
   public:
      explicit HZZFourLElectronSelector(const edm::ParameterSet&);
      ~HZZFourLElectronSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      //virtual void beginStream(edm::StreamID) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) /*override*/;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) /*override*/;
      bool leptonCrosscleaned(const pat::Electron& ele,const edm::Handle<pat::MuonCollection>& muVec);
      //bool passBDT(const double fSCeta, const double pt, const double BDT);
      //virtual void endStream() override;
      
      std::vector<pat::Electron>* looselist_;
      std::vector<pat::Electron>* looseSIPlist_;
      std::vector<pat::Electron>* tightlist_;

      int verbosity_;

      const edm::InputTag vertexTag_;
      const edm::InputTag electronTag_;
      const edm::InputTag tightSIPMuonTag_;


      const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      const edm::EDGetTokenT<pat::MuonCollection> tightSIPMuonToken_;

      const std::string looseEleOutputColl_;
      const std::string looseSIPEleOutputColl_;
      const std::string tightEleOutputColl_;
};
#endif
