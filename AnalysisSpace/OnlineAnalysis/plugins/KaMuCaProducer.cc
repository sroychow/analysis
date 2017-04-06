#include <iostream>
#include <algorithm>

#include "AnalysisSpace/OnlineAnalysis/plugins/KaMuCaProducer.h"


//
// constructors and destructor
//
KaMuCaProducer::KaMuCaProducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  isSync_(iConfig.getParameter<bool>("isSync")),
  muonToken_(consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc")))
{
  if (isMC_) kaMuCa_ = new KalmanMuonCalibrator("MC_80X_13TeV");
  else kaMuCa_ = new KalmanMuonCalibrator("DATA_80X_13TeV");
  produces<std::vector<pat::Muon> >(); 
}


KaMuCaProducer::~KaMuCaProducer()
{
}


// ------------ method called to produce the data  ------------
void
KaMuCaProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  std::vector<pat::Muon>*   calibmuons = new std::vector<pat::Muon>();
  for( const auto& mu : *muons ) {
    double muPt = mu.pt();
    double corrPt =  muPt;
    double pterr = mu.muonBestTrack()->ptError();
    double corrPtError = pterr;
    if(mu.muonBestTrackType() != 1 || mu.pt() >= 200.0) continue;
    if(isMC_) {
      corrPt = kaMuCa_->getCorrectedPt(muPt, mu.eta(), mu.phi(), mu.charge());
      double smearedPt;
      if(isSync_) smearedPt = kaMuCa_->smearForSync(corrPt, mu.eta());
      else smearedPt = kaMuCa_->smear(corrPt, mu.eta());
      corrPtError = smearedPt * kaMuCa_->getCorrectedError(smearedPt, mu.eta(), pterr/smearedPt );//pr error after smearing
      corrPt = smearedPt;
    } else {
      corrPt = kaMuCa_->getCorrectedPt(muPt, mu.eta(), mu.phi(), mu.charge());
      corrPtError = corrPt * kaMuCa_->getCorrectedError(corrPt, mu.eta(), pterr/corrPt );
    } 
    pat::Muon calibmu = mu;
    calibmu.addUserFloat("correctedPtError",corrPtError);
    TLorentzVector tempP4;
    tempP4.SetPtEtaPhiM(corrPt, mu.eta(), mu.phi(), mu.mass());
    calibmu.setP4(reco::Particle::PolarLorentzVector(tempP4.Pt(),tempP4.Eta(),tempP4.Phi(),tempP4.M()));
    calibmuons->push_back(calibmu); 
  }
  // add the muons to the event output
  std::auto_ptr<std::vector<pat::Muon> > ptr(calibmuons);
  iEvent.put(ptr);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
KaMuCaProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
KaMuCaProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  KaMuCaProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  KaMuCaProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
KaMuCaProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
KaMuCaProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KaMuCaProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KaMuCaProducer);
