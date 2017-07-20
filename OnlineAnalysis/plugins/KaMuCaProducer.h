#ifndef __AnalysisSpace_OnlineAnalysis_KaMuCaProducer_h
#define __AnalysisSpace_OnlineAnalysis_KaMuCaProducer_h
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

//Root
#include "TLorentzVector.h"

// Kalman Muon Corrections 
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

//
// class declaration
//

class KaMuCaProducer : public edm::stream::EDProducer<> {
   public:
      explicit KaMuCaProducer(const edm::ParameterSet&);
      ~KaMuCaProducer();

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
      bool isMC_;
      bool isSync_;
      const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      //Kalman Muon Calibration 
      KalmanMuonCalibrator *kaMuCa_;
};
#endif
