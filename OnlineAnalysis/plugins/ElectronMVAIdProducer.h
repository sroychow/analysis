#ifndef __AnalysisSpace_OnlineAnalysis_ElectronMVAIdProducer_h
#define __AnalysisSpace_OnlineAnalysis_ElectronMVAIdProducer_h
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

//Root
#include "TLorentzVector.h"

// Kalman Muon Corrections 
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

//
// class declaration
//

class ElectronMVAIdProducer : public edm::stream::EDProducer<> {
   public:
      explicit ElectronMVAIdProducer(const edm::ParameterSet&);
      ~ElectronMVAIdProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      bool verbosity_;
      const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      const edm::EDGetTokenT<edm::View<reco::GsfElectron> > gsfelectronTokenMVAId_;
      // ID decisions objects
      //edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      //edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  
      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      //edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
};
#endif
