#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TVector2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "AnalysisSpace/OnlineAnalysis/plugins/ZCandidateProducer.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"


//
// constructors and destructor
//
ZCandidateProducer::ZCandidateProducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  tightIsoEleFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightIsoElectronFSRSrc")),
  tightIsoMuFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightIsoMuonFSRSrc")),
  tightIsoElectronFSRToken_(consumes<std::vector<vhtm::IsoElectron>>(tightIsoEleFSRpairTag_)),
  tightIsoMuonFSRToken_(consumes<std::vector<vhtm::IsoMuon>>(tightIsoMuFSRpairTag_)),
  ZeeColl_(iConfig.getUntrackedParameter<std::string>("ZEEColl", "ZToeeList")),
  ZmumuColl_(iConfig.getUntrackedParameter<std::string>("ZMuMuColl", "ZTomumuList"))
{
  produces<std::vector<vhtm::Zee>>(ZeeColl_);
  produces<std::vector<vhtm::Zmumu>>(ZmumuColl_);
}


ZCandidateProducer::~ZCandidateProducer()
{
}


// ------------ method called to produce the data  ------------
void
ZCandidateProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ZeeVec_->clear();
  ZmumuVec_->clear();
  
  //using namespace edm;
  
  edm::Handle<std::vector<vhtm::IsoElectron>> tightIsoEle;
  iEvent.getByToken(tightIsoElectronFSRToken_, tightIsoEle);

  edm::Handle<std::vector<vhtm::IsoMuon>> tightIsoMu;
  iEvent.getByToken(tightIsoMuonFSRToken_, tightIsoMu);
  
  ZeeSelector(*tightIsoEle,*ZeeVec_);
  ZmumuSelector(*tightIsoMu, *ZmumuVec_);

  std::auto_ptr<std::vector<vhtm::Zee>> pv1(new std::vector<vhtm::Zee>(*ZeeVec_));
  iEvent.put(pv1,ZeeColl_);
  std::auto_ptr<std::vector<vhtm::Zmumu>> pv2(new std::vector<vhtm::Zmumu>(*ZmumuVec_));
  iEvent.put(pv2,ZmumuColl_);
}

void 
ZCandidateProducer::ZmumuSelector(const std::vector<vhtm::IsoMuon>& lepPhotonPairVec, std::vector<vhtm::Zmumu>& candList)
{
  //std::cout << "Zmm" << std::endl;
  for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
    const auto& ip = lepPhotonPairVec[i];
    const TLorentzVector& lep1P4 = HZZ4lUtil::getP4(ip.mu);
    TLorentzVector lep1fsrP4;
    lep1fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
    if(!ip.hasfsr) lep1fsrP4 = HZZ4lUtil::getP4(ip.fsr);    
    for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
      const auto& jp = lepPhotonPairVec[j];

      if (ip.mu.charge() + jp.mu.charge() != 0) continue; 

      const TLorentzVector& lep2P4 = HZZ4lUtil::getP4(jp.mu);
      TLorentzVector lep2fsrP4;
      lep2fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
      if(!jp.hasfsr) lep2fsrP4 = HZZ4lUtil::getP4(jp.fsr);

      double zM = (lep1P4 + lep2P4 + lep1fsrP4 + lep2fsrP4).M();
      //if(zM <= 12. && zM >= 120.)     continue;
      vhtm::Zmumu Ztemp;
      Ztemp.lep1 = ip.mu;
      Ztemp.lep2 = jp.mu;
      Ztemp.fsrl1 = ip.fsr;
      Ztemp.fsrl2 = jp.fsr;
      Ztemp.lep1hasfsr = ip.hasfsr;
      Ztemp.lep2hasfsr = jp.hasfsr;
      Ztemp.lep1Iso = ip.relIso;
      Ztemp.lep2Iso = jp.relIso;
      Ztemp.mass = zM;
      Ztemp.massdiff = std::fabs(HZZ4lUtil::MZnominal - zM);
      candList.push_back(Ztemp);
    }
  }
}
void 
ZCandidateProducer::ZeeSelector(const std::vector<vhtm::IsoElectron>& lepPhotonPairVec, std::vector<vhtm::Zee>& candList)
{
  //std::cout << "Zee" << std::endl;
  for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
    const auto& ip = lepPhotonPairVec[i];
    const TLorentzVector& lep1P4 = HZZ4lUtil::getP4(ip.ele);
    TLorentzVector lep1fsrP4;
    lep1fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
    if(!ip.hasfsr) lep1fsrP4 = HZZ4lUtil::getP4(ip.fsr);    
    for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
      const auto& jp = lepPhotonPairVec[j];

      if (ip.ele.charge() + jp.ele.charge() != 0) continue; 

      const TLorentzVector& lep2P4 = HZZ4lUtil::getP4(jp.ele);
      TLorentzVector lep2fsrP4;
      lep2fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
      if(!jp.hasfsr) lep2fsrP4 = HZZ4lUtil::getP4(jp.fsr);

      double zM = (lep1P4 + lep2P4 + lep1fsrP4 + lep2fsrP4).M();
      if(zM <= 12. && zM >= 120.)     continue;
      vhtm::Zee Ztemp;
      Ztemp.lep1 = ip.ele;
      Ztemp.lep2 = jp.ele;
      Ztemp.fsrl1 = ip.fsr;
      Ztemp.fsrl2 = jp.fsr;
      Ztemp.lep1hasfsr = ip.hasfsr;
      Ztemp.lep2hasfsr = jp.hasfsr;
      Ztemp.lep1Iso = ip.relIso;
      Ztemp.lep2Iso = jp.relIso;
      Ztemp.mass = zM;
      Ztemp.massdiff = std::fabs(HZZ4lUtil::MZnominal - zM);
      candList.push_back(Ztemp);
    }
  }
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ZCandidateProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ZCandidateProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  ZCandidateProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  ZCandidateProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
ZCandidateProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  ZeeVec_= new std::vector<vhtm::Zee>();
  ZmumuVec_= new std::vector<vhtm::Zmumu>();
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
ZCandidateProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZCandidateProducer);
