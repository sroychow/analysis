#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVector3.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "AnalysisSpace/OnlineAnalysis/plugins/HZZFourLElectronSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"




HZZFourLElectronSelector::HZZFourLElectronSelector(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc", edm::InputTag("selectedPatElectrons"))),
  tightSIPMuonTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightSIPMuonSrc",edm::InputTag("tightSIPMuonVector"))),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  electronToken_(consumes<pat::ElectronCollection>(electronTag_)),
  tightSIPMuonToken_(consumes<pat::MuonCollection>(tightSIPMuonTag_)),
  looseEleOutputColl_(iConfig.getUntrackedParameter<std::string>("looseElecoll","looseElectronVector")),
  looseSIPEleOutputColl_(iConfig.getUntrackedParameter<std::string>("looseSIPElecoll","looseSIPElectronVector")),
  tightEleOutputColl_(iConfig.getUntrackedParameter<std::string>("tightElecoll","tightElectronVector"))
{
  produces<std::vector<pat::Electron>>(looseEleOutputColl_);
  produces<std::vector<pat::Electron>>(looseSIPEleOutputColl_);
  produces<std::vector<pat::Electron>>(tightEleOutputColl_);
}


HZZFourLElectronSelector::~HZZFourLElectronSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HZZFourLElectronSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  looselist_ = new std::vector<pat::Electron>();
  looseSIPlist_ = new std::vector<pat::Electron>();
  tightlist_ = new std::vector<pat::Electron>();
}
//
// member functions
//

// ------------ method called to produce the data  ------------
void
HZZFourLElectronSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //std::cout << "Entering HZZFourLElectronSelector::produce" << std::endl;
  // using namespace edm;
   looselist_->clear();
   looseSIPlist_->clear();
   tightlist_->clear();

   edm::Handle<pat::ElectronCollection> electrons;
   bool found = iEvent.getByToken(electronToken_, electrons);
   edm::LogInfo("HZZFourLElectronSelector") << "Total # of Electrons: " << electrons->size();
   
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(tightSIPMuonToken_, muons);
   //std::cout << "Total # of Electrons: " <<electrons->size()<<std::endl;

   if (found && electrons.isValid()) {

     edm::Handle<reco::VertexCollection> primaryVertices;
     iEvent.getByToken(vertexToken_, primaryVertices);

     edm::LogInfo("HZZFourLElectronSelector") << "Total # PAT Electrons: " << electrons->size();

     /*
     //Loose Electrons: pT > 7, |eta| < 2.5, dxy< 0.5, dz < 1, 
     //Tight Electrons: Loose Electrons + non triggering MVA ID
     //In CMSSW_7_6_X , we don't recompute the ID and we directly retrieve the variable from the 
     //pat::Electron as userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")
     */
     for (const auto& ele: *electrons) {

       bool hasGsfTrack = ele.gsfTrack().isNonnull() ? true : false;
       if (ele.pt() <= 7 ) continue;
       if (fabs(ele.eta()) >= 2.5) continue;
       
       double dxyWrtPV = -99.;
       double dzWrtPV = -99.;
       if (hasGsfTrack) {
	 reco::GsfTrackRef tk = ele.gsfTrack();
	 if (primaryVertices.isValid()) {
	   const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
	   dxyWrtPV = tk->dxy(vit.position());
	   dzWrtPV  = tk->dz(vit.position());
	 }
       }
       if(dxyWrtPV >= 0.5) continue;
       if(dzWrtPV >= 1.) continue;
       looselist_->push_back(ele);
       if(std::fabs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D)) >= 4.)    continue;
       if(!leptonCrosscleaned(ele,muons))   continue;
       looseSIPlist_->push_back(ele);
       if(HZZ4lUtil::passBDT(std::fabs(ele.superCluster()->eta()), 
                  ele.pt(), 
                  ele.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")))
         tightlist_->push_back(ele);
     }
     std::auto_ptr<std::vector<pat::Electron>> pv1(new std::vector<pat::Electron>(*looselist_));       
     iEvent.put(pv1,looseEleOutputColl_);
     std::auto_ptr<std::vector<pat::Electron>> pv2(new std::vector<pat::Electron>(*looseSIPlist_));       
     iEvent.put(pv2,looseSIPEleOutputColl_);
     std::auto_ptr<std::vector<pat::Electron>> pv3(new std::vector<pat::Electron>(*tightlist_));       
     iEvent.put(pv3,tightEleOutputColl_);
   }
  //std::cout << "Leaving HZZFourLElectronSelector::produce" << std::endl;
}

bool HZZFourLElectronSelector::leptonCrosscleaned(const pat::Electron& ele,const edm::Handle<pat::MuonCollection>& muVec) {
  bool flag = true;
  std::cout << "Entering Cross clean" << std::endl;
  for (const auto& mu: *muVec) {
    if (HZZ4lUtil::getP4(ele).DeltaR(HZZ4lUtil::getP4(mu)) < 0.05) {
      flag = false;
      break;
    }
  }
  return flag;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
/*void
HZZFourLElectronSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HZZFourLElectronSelector::endStream() {
}
*/
// ------------ method called when starting to processes a run  ------------
/*
void
HZZFourLElectronSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HZZFourLElectronSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------


 
// ------------ method called when ending the processing of a luminosity block  ------------

void
HZZFourLElectronSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HZZFourLElectronSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZFourLElectronSelector);
