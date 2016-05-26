// -*- C++ -*-
//
// Package:    AnalysisSpace/TriggerFilter
// Class:      TriggerFilter
// 
/**\class TriggerFilter TriggerFilter.cc AnalysisSpace/TriggerFilter/plugins/TriggerFilter.cc

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
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "TPRegexp.h"
//
// class declaration
//

class TriggerFilter : public edm::stream::EDFilter<> {
   public:
      explicit TriggerFilter(const edm::ParameterSet&);
      ~TriggerFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override;
      virtual void endRun(edm::Run const& , edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::InputTag hltTag_;
      const std::vector<std::string> hltPathsOfInterest_;
      std::vector<std::string> matchedPathList_;
      HLTConfigProvider hltConfig_;
      const edm::EDGetTokenT<edm::TriggerResults> hltToken_;
      HLTPrescaleProvider hltPrescaleProvider_;
      int nTrigpassEvt;
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
TriggerFilter::TriggerFilter(const edm::ParameterSet& iConfig) :
  hltTag_(iConfig.getUntrackedParameter<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults","","HLT"))),
  hltPathsOfInterest_(iConfig.getParameter<std::vector<std::string> >("hltPathsOfInterest")),
  hltToken_(consumes<edm::TriggerResults>(hltTag_)),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
   //now do what ever initialization is needed
   nTrigpassEvt = 0;

}


TriggerFilter::~TriggerFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool isTriggered = false;
  edm::Handle<edm::TriggerResults> triggerResults;
  bool found = iEvent.getByToken(hltToken_, triggerResults);
  if (found && triggerResults.isValid()) {
    for (auto path: matchedPathList_) {
      int fired = -1;
      unsigned int index = hltConfig_.triggerIndex(path);
      if (index < triggerResults->size()) {
        fired = (triggerResults->accept(index)) ? 1 : 0;
      }
      else {
	edm::LogInfo("TriggerFilter") << "Requested HLT path \"" << path << "\" does not exist";
      }
      int prescale = -1;
      const int prescaleSet = hltPrescaleProvider_.prescaleSet(iEvent, iSetup);
      //if (hltConfig_.prescaleSet(iEvent, iSetup) < 0) {
      if ( prescaleSet < 0 ) {
	edm::LogError("TriggerBlock") << "The prescale set index number could not be obtained for HLT path: "
                                      << path;
      }
      else {
        prescale = hltConfig_.prescaleValue(prescaleSet, path);
      }
      if(fired == 1 && prescale == 1)    {
        isTriggered = true;
        nTrigpassEvt++;
        break;
      }
    }
  }
  return isTriggered;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TriggerFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TriggerFilter::endStream() {
  std::cout << "Events Passing Trigger(endStream) = " << nTrigpassEvt << std::endl;
}

// ------------ method called when starting to processes a run  ------------

void
TriggerFilter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{ 
  bool changed = true;
  if (hltConfig_.init(iRun, iSetup, hltTag_.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerBlock") << "HLT config with process name "
				 << hltTag_.process() << " successfully extracted";
    matchedPathList_.clear();
    const std::vector<std::string>& pathList = hltConfig_.triggerNames();
    for (const std::string& path: pathList) {
      if (hltPathsOfInterest_.size()) {
        int nmatch = 0;
        for (const std::string& kt: hltPathsOfInterest_)
          nmatch += TPRegexp(kt).Match(path);
        if (!nmatch) continue;
      }
      matchedPathList_.push_back(path);
    }
  }
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerBlock") << "Error! HLT config extraction with process name "
				  << hltTag_.process() << " failed";
    // In this case, all access methods will return empty values!
  }
}

// ------------ method called when ending the processing of a run  ------------

void
TriggerFilter::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  std::cout << "Events Passing Trigger(endRun) = " << nTrigpassEvt << std::endl;
}

 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TriggerFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TriggerFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TriggerFilter);
