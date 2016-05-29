// -*- C++ -*-
#include <iostream>
#include<iomanip>
#include "TTree.h"
#include "TClass.h"
using std::setw;

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "AnalysisSpace/OnlineAnalysis/plugins/HZZFourLJetSelector.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

//
// Package:    HZZFourLJetSelector/MuonBlock
// Class:      HZZFourLJetSelector
// 
/**\class HZZFourLJetSelector MuonBlock.cc MuonBlock/MuonBlock/plugins/MuonBlock.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

//
// static data member definitions
//

//
// constructors and destructor
//

PFJetIDSelectionFunctor pfjetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
PFJetIDSelectionFunctor pfjetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

HZZFourLJetSelector::HZZFourLJetSelector(const edm::ParameterSet& iConfig):
  tightIsoEleFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightIsoElectronFSRSrc")),
  tightIsoMuFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("tightIsoMuonFSRSrc")),
  jetTag_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc", edm::InputTag("selectedPatJets"))), 
  tightIsoElectronFSRToken_(consumes<std::vector<vhtm::IsoElectron>>(tightIsoEleFSRpairTag_)),
  tightIsoMuonFSRToken_(consumes<std::vector<vhtm::IsoMuon>>(tightIsoMuFSRpairTag_)),
  jetToken_(consumes<pat::JetCollection>(jetTag_)),
  looseJetOutputColl_(iConfig.getUntrackedParameter<std::string>("looseJetColl","selectedlooseJets"))
{
  produces<std::vector<pat::Jet>>(looseJetOutputColl_);
}


HZZFourLJetSelector::~HZZFourLJetSelector()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HZZFourLJetSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//using namespace edm;
	looseJetVec_->clear();


        edm::Handle<std::vector<vhtm::IsoElectron>> tightIsoEle;
 	iEvent.getByToken(tightIsoElectronFSRToken_, tightIsoEle);

        edm::Handle<std::vector<vhtm::IsoMuon>> tightIsoMu;
        iEvent.getByToken(tightIsoMuonFSRToken_, tightIsoMu);	
    
        edm::Handle<pat::JetCollection> jetColl;
	bool jetFound = iEvent.getByToken(jetToken_, jetColl);
        
        if(jetFound) {
          std::vector<TLorentzVector>   p4vector;
          for(auto& e: *tightIsoEle) {
            p4vector.push_back(HZZ4lUtil::getP4(e.ele));
            if(e.hasfsr)    p4vector.push_back(HZZ4lUtil::getP4(e.fsr));
          }
          for(auto& m: *tightIsoMu) {
            p4vector.push_back(HZZ4lUtil::getP4(m.mu));
            if(m.hasfsr)    p4vector.push_back(HZZ4lUtil::getP4(m.fsr));
          }
          for(auto& j: *jetColl) {
	    if(j.pt() <= 30)       continue;
            if(std::fabs(j.eta()) >= 4.7)   continue;
            if(!jetLeptonCleaning(HZZ4lUtil::getP4(j), p4vector,0.4))      continue;
            if (!isLooseJet(j))    continue;
            looseJetVec_->push_back(j);
          }
	
        }
	std::auto_ptr<std::vector<pat::Jet>> pv1(new std::vector<pat::Jet>(*looseJetVec_));       
	iEvent.put(pv1,looseJetOutputColl_);
}

bool HZZFourLJetSelector::jetLeptonCleaning(const TLorentzVector& jetP4, const std::vector<TLorentzVector>& lepfsrp4vec,double dR) 
{
  for(auto& lp4: lepfsrp4vec) {
    if(lp4.DeltaR(jetP4) <= 0.4)      return false;
  }
  return true;
}

bool HZZFourLJetSelector::isLooseJet(const pat::Jet& jet) {
  bool centralCut = (std::fabs(jet.eta()) <= 2.4) 
      ? (jet.chargedHadronEnergyFraction() > 0 && 
  	 jet.chargedMultiplicity() > 0 && 
  	 jet.chargedEmEnergyFraction() < 0.99)
      : true;
    
  return (jet.neutralHadronEnergyFraction() < 0.99 && 
  	  jet.neutralEmEnergyFraction() < 0.99 &&
  	  (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1 &&
  	  jet.muonEnergyFraction() < 0.8 &&
  	  centralCut);
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HZZFourLJetSelector::beginJob(/*edm::StreamID*/)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HZZFourLJetSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
HZZFourLJetSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HZZFourLJetSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------

void
HZZFourLJetSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  looseJetVec_ = new std::vector<pat::Jet>();
}
 
// ------------ method called when ending the processing of a luminosity block  ------------

void
HZZFourLJetSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HZZFourLJetSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZFourLJetSelector);
