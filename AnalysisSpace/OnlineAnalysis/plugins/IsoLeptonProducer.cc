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

#include "AnalysisSpace/OnlineAnalysis/plugins/IsoLeptonProducer.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"


//
// constructors and destructor
//
IsoLeptonProducer::IsoLeptonProducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  fsrTag_(iConfig.getUntrackedParameter<edm::InputTag>("fsrCol",edm::InputTag("packedPFCandidates"))), 
  looseSIPEleFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("looseSIPElectronFSRSrc")),
  looseSIPMuFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("looseSIPMuonFSRSrc")),
  fixedGridRhoFastjetAllTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllTag",edm::InputTag("fixedGridRhoFastjetAll"))),
  fsrToken_(consumes<pat::PackedCandidateCollection>(fsrTag_)),
  looseElectronFSRToken_(consumes<EleFSRCollection>(looseSIPEleFSRpairTag_)),
  looseMuonFSRToken_(consumes<MuFSRCollection>(looseSIPMuFSRpairTag_)),
  fixedGridRhoFastjetAllToken_(consumes<double>(fixedGridRhoFastjetAllTag_)),
  tightSIPEleFSRColl_(iConfig.getUntrackedParameter<std::string>("tightSIPEleFSRColl", "tightEleFSRPairList")),
  tightSIPMuFSRColl_(iConfig.getUntrackedParameter<std::string>("tightSIPMuFSRColl", "tightMuFSRPairList"))
{
  produces<std::vector<vhtm::IsoElectron>>(tightSIPEleFSRColl_);
  produces<std::vector<vhtm::IsoMuon>>(tightSIPMuFSRColl_);
}


IsoLeptonProducer::~IsoLeptonProducer()
{
}


// ------------ method called to produce the data  ------------
void
IsoLeptonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  tightEleFSRpairVec_->clear();
  tightMuFSRpairVec_->clear();
  
  //using namespace edm;
  
  edm::Handle<pat::PackedCandidateCollection> fsr;
  iEvent.getByToken(fsrToken_,fsr);
  edm::Handle<EleFSRCollection> looseSIPEleFSRpair;
  iEvent.getByToken(looseElectronFSRToken_,looseSIPEleFSRpair);
  edm::Handle<MuFSRCollection> looseSIPMuFSRpair;
  iEvent.getByToken(looseMuonFSRToken_,looseSIPMuFSRpair);
  
  edm::Handle<double> fixedGridRhoFastjetAll;
  iEvent.getByToken(fixedGridRhoFastjetAllToken_,fixedGridRhoFastjetAll);
  const double evrho = *fixedGridRhoFastjetAll;
  
  for(auto& e : *looseSIPEleFSRpair) {
     //tight electron cut
     if(!HZZ4lUtil::passBDT(std::fabs(e.first.superCluster()->eta()), e.first.pt(), 
                  e.first.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")))  continue;
     double iso = HZZ4lUtil::computeElectronReliso(e.first, fsr,evrho);
     if( iso >= 0.35)      continue;
     vhtm::IsoElectron ietemp;
     ietemp.ele = e.first;
     if(!e.second.empty())  {
       ietemp.hasfsr = true;//default is false
       ietemp.fsr = e.second.at(0);
     }
     ietemp.relIso = iso;
     tightEleFSRpairVec_->push_back(ietemp);
  }

  for(auto& m : *looseSIPMuFSRpair) {
     //tight mu cut
     if(!m.first.isPFMuon())    continue; 
     double iso = HZZ4lUtil::computeMuonReliso(m.first, fsr);
     if( iso >= 0.35)      continue;
     vhtm::IsoMuon imtemp;
     imtemp.mu = m.first;
     if(!m.second.empty())  {
       imtemp.hasfsr = true;//default is false
       imtemp.fsr = m.second.at(0);
     }
     imtemp.relIso = iso;
     tightMuFSRpairVec_->push_back(imtemp);
  }
  

  std::auto_ptr<std::vector<vhtm::IsoElectron>> pv1(new std::vector<vhtm::IsoElectron>(*tightEleFSRpairVec_));
  iEvent.put(pv1,tightSIPEleFSRColl_);
  std::auto_ptr<std::vector<vhtm::IsoMuon>> pv2(new std::vector<vhtm::IsoMuon>(*tightMuFSRpairVec_));
  iEvent.put(pv2,tightSIPMuFSRColl_);
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
IsoLeptonProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
IsoLeptonProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  IsoLeptonProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  IsoLeptonProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
IsoLeptonProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  tightEleFSRpairVec_= new std::vector<vhtm::IsoElectron>();
  tightMuFSRpairVec_ = new std::vector<vhtm::IsoMuon>();
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
IsoLeptonProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
IsoLeptonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoLeptonProducer);
