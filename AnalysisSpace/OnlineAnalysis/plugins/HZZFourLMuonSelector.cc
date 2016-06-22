// -*- C++ -*-
#include <iostream>
#include<iomanip>
#include "TTree.h"
#include "TClass.h"
using std::setw;

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "AnalysisSpace/OnlineAnalysis/plugins/HZZFourLMuonSelector.h"

#include "AnalysisSpace/OnlineAnalysis/plugins/HZZFourLMuonSelector.h"
//
// Package:    HZZFourLMuonSelector/HZZFourLMuonSelector
// Class:      HZZFourLMuonSelector
// 
/**\class HZZFourLMuonSelector HZZFourLMuonSelector.cc HZZFourLMuonSelector/HZZFourLMuonSelector/plugins/HZZFourLMuonSelector.cc

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
HZZFourLMuonSelector::HZZFourLMuonSelector(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc", edm::InputTag("selectedPatMuons"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  muonToken_(consumes<pat::MuonCollection>(muonTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  looseMuOutputColl_(iConfig.getUntrackedParameter<std::string>("looseMucoll","looseMuonVector")),
  looseSIPMuOutputColl_(iConfig.getUntrackedParameter<std::string>("looseSIPMucoll","looseSIPMuonVector")),
  tightMuOutputColl_(iConfig.getUntrackedParameter<std::string>("tightMucoll","tightMuonVector"))
{
  produces<std::vector<pat::Muon>>(looseMuOutputColl_);
  produces<std::vector<pat::Muon>>(looseSIPMuOutputColl_);
  produces<std::vector<pat::Muon>>(tightMuOutputColl_);
}


HZZFourLMuonSelector::~HZZFourLMuonSelector()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HZZFourLMuonSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//using namespace edm;
  //std::cout << "Entering HZZFourLMuonSelector::produce" << std::endl;
	looselist_->clear();
        looseSIPlist_->clear();
	tightlist_->clear();
	fnMuon_ = 0;

	edm::Handle<pat::MuonCollection> muons;
	bool found = iEvent.getByToken(muonToken_, muons);


	if (found ) {
	  edm::Handle<reco::VertexCollection> primaryVertices;
	  iEvent.getByToken(vertexToken_, primaryVertices);


	  edm::LogInfo("HZZFourLMuonSelector") << "Total # of Muons: " << muons->size();
	  /*
          Loose Muons: pT > 5, 
                       |eta| < 2.4, 
                       dxy< 0.5, 
                       dz < 1, 
                       (isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && 
                       muonBestTrackType!=2
                       dxy and dz are defined w.r.t. the PV and using the muonBestTrack, 
                       e.g, dxy = fabs(l.muonBestTrack()->dxy(PV->position()))
                       **Note that non-global tracker muons must be arbitrated (numberOfMatches>0) and 
                       that muons with muonBestTrackType == 2 (standalone) are discarded even if they 
                       are marked as global or traker muons. 
          Tight Muons: as Loose Muons+ PF Muon 
          */
	  for( const auto& mu : *muons ) {
            if(mu.pt() <= 5.) continue;
            if(std::fabs(mu.eta()) > 2.4) continue;
	    reco::TrackRef tk = mu.muonBestTrack();
	    double dxyWrtPV = -99.;
	    double dzWrtPV = -99.;
            bool highPtid = false;
	    if (primaryVertices.isValid()) {
	      edm::LogInfo("HZZFourLMuonSelector") << "Total # Primary Vertices: " << primaryVertices->size();

	      const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
	      dxyWrtPV = tk->dxy(vit.position());
	      dzWrtPV  = tk->dz(vit.position());
              highPtid = isTrackerHighPt(mu,vit) && mu.pt() > 200;
	    }
	    if(std::fabs(dxyWrtPV) >= 0.5 )      continue;
            if(std::fabs(dzWrtPV) >= 1.)         continue;
	    bool quality = (mu.isGlobalMuon()
                           || ( mu.isTrackerMuon() && mu.numberOfMatches() >= 0)) 
                           && mu.muonBestTrackType()!=2 ;
            if(!quality) continue;
            looselist_->push_back(mu);
            if(std::fabs(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D)) >= 4.)    continue;
            looseSIPlist_->push_back(mu);
            bool isTight = mu.isPFMuon() || highPtid;
            if(!isTight)    continue; 
            tightlist_->push_back(mu);
	  }
	}
        //std::cout<<"No. of loose muons"<<looselist_->size()<<std::endl;
	std::auto_ptr<std::vector<pat::Muon>> pv1(new std::vector<pat::Muon>(*looselist_));       
	iEvent.put(pv1,looseMuOutputColl_);
	std::auto_ptr<std::vector<pat::Muon>> pv2(new std::vector<pat::Muon>(*looseSIPlist_));       
	iEvent.put(pv2,looseSIPMuOutputColl_);
	std::auto_ptr<std::vector<pat::Muon>> pv3(new std::vector<pat::Muon>(*tightlist_));       
	iEvent.put(pv3,tightMuOutputColl_);
  //std::cout << "Leaving HZZFourLMuonSelector::produce" << std::endl;
}

bool HZZFourLMuonSelector::isTrackerHighPt(const pat::Muon & mu, const reco::Vertex & primaryVertex) {
  //reco::TrackRef
  const auto& bestrkRef = mu.muonBestTrack();
  const auto& intrkRef = mu.innerTrack();
  if(!bestrkRef.isNonnull() || !intrkRef.isNonnull() )      return false;		
  return ( mu.numberOfMatchedStations() > 1 
                         && (bestrkRef->ptError()/bestrkRef->pt()) < 0.3 
                         && std::abs(bestrkRef->dxy(primaryVertex.position())) < 0.2 
                         && std::abs(bestrkRef->dz(primaryVertex.position())) < 0.5 
                         && intrkRef->hitPattern().numberOfValidPixelHits() > 0 
                         && intrkRef->hitPattern().trackerLayersWithMeasurement() > 5 );

}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HZZFourLMuonSelector::beginJob(/*edm::StreamID*/)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HZZFourLMuonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
HZZFourLMuonSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HZZFourLMuonSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------

void
HZZFourLMuonSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  looselist_ = new std::vector<pat::Muon>();
  looseSIPlist_ = new std::vector<pat::Muon>();
  tightlist_ = new std::vector<pat::Muon>();
}
 
// ------------ method called when ending the processing of a luminosity block  ------------

void
HZZFourLMuonSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HZZFourLMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZFourLMuonSelector);
