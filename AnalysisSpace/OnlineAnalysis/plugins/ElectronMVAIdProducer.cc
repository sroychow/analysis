#include <iostream>
#include <algorithm>

#include "AnalysisSpace/OnlineAnalysis/plugins/ElectronMVAIdProducer.h"


//
// constructors and destructor
//
ElectronMVAIdProducer::ElectronMVAIdProducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  gsfelectronTokenMVAId_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))) 
{
  produces<std::vector<pat::Electron> >(); 
}


ElectronMVAIdProducer::~ElectronMVAIdProducer()
{
}


// ------------ method called to produce the data  ------------
void
ElectronMVAIdProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle<edm::View<reco::GsfElectron> > gsfelectrons;
  iEvent.getByToken(gsfelectronTokenMVAId_, gsfelectrons);

  //edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  //edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
  //iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  //iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > mvaValues;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  //edm::Handle<edm::ValueMap<int> > mvaCategories;
  //iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

  //output electrons
  std::vector<pat::Electron>* outelectrons = new std::vector<pat::Electron>();
  //auto gsfit = gsfelectrons->begin();
  unsigned int gsfeleidx = 0;
  for (const pat::Electron& v: *electrons) {
      //storing of ele id decisions
      const auto gsfel = gsfelectrons->ptrAt(gsfeleidx);
      //electron.passMediumId = (*medium_id_decisions)[gsfelgsfel];
      //electron.passTightId = (*tight_id_decisions)[gsfel];
      //double eleBDT = (*mvaValues)[gsfel]; 
      //electron.mvaCategory = (*mvaCategories)[gsfel];
      gsfeleidx++;
      pat::Electron oele = v;
      oele.addUserFloat("electronBDT",(*mvaValues)[gsfel]);
    outelectrons->push_back(oele); 
  }
  // add the muons to the event output
  std::auto_ptr<std::vector<pat::Electron> > ptr(outelectrons);
  iEvent.put(ptr);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ElectronMVAIdProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ElectronMVAIdProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  ElectronMVAIdProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  ElectronMVAIdProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
ElectronMVAIdProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
ElectronMVAIdProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronMVAIdProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMVAIdProducer);
