// -*- C++ -*-
//
// Package:    AnalysisSpace/ZZFilter
// Class:      ZZFilter
// 
/**\class ZZFilter ZZFilter.cc AnalysisSpace/ZZFilter/plugins/ZZFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Suvankar Roy Chowdhury
//         Created:  Sat, 07 May 2016 10:30:39 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"

//
// class declaration
//

class ZZFilter : public edm::stream::EDFilter<> {
   public:
      explicit ZZFilter(const edm::ParameterSet&);
      ~ZZFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::InputTag ZZCandVecTag_;
      const edm::EDGetTokenT<std::vector<vhtm::ZZcandidate>> ZZVecToken_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZZFilter::ZZFilter(const edm::ParameterSet& iConfig) :
  ZZCandVecTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZZCandVec")),
  ZZVecToken_(consumes<std::vector<vhtm::ZZcandidate>>(ZZCandVecTag_))
{
   //now do what ever initialization is needed

}


ZZFilter::~ZZFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZZFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<std::vector<vhtm::ZZcandidate>> ZZvector_;
   iEvent.getByToken(ZZVecToken_, ZZvector_);
   if(ZZvector_->empty())    return false;
   return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ZZFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ZZFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ZZFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ZZFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ZZFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ZZFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZZFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZZFilter);
