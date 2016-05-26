#ifndef __AnalysisSpace_OnlineAnalysis_ZZCandidateProducer_h
#define __AnalysisSpace_OnlineAnalysis_ZZCandidateProducer_h
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

class ZZCandidateProducer : public edm::stream::EDProducer<> {
   public:
      explicit ZZCandidateProducer(const edm::ParameterSet&);
      ~ZZCandidateProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      
      template <typename T1,typename T2>
       int ZZSelector(const T1& Zcand1, const T2& Zcand2, int& whichZ1cand,double& m4l, bool sameFlavour = false);

      int verbosity_;
      const edm::InputTag ZeeTag_;
      const edm::InputTag ZmumuTag_;
      

      const edm::EDGetTokenT<std::vector<vhtm::Zee>> ZeeToken_;
      const edm::EDGetTokenT<std::vector<vhtm::Zmumu>> ZmumuToken_;

      const std::string ZZColl_;

      std::vector<vhtm::ZZcandidate>* ZZVec_;
};
#endif
