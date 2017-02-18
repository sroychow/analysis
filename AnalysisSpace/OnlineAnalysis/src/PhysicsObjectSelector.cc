#include <iostream>
#include <algorithm>

#include "AnalysisSpace/OnlineAnalysis/src/PhysicsObjectSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/LeptonIsoCalculator.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

PhysicsObjectSelector::PhysicsObjectSelector(const bool verbose) {
  verbosity_ = verbose;
  looseMulist_ = new std::vector<pat::Muon>();
  looseMuSIPlist_ = new std::vector<pat::Muon>();
  tightMulist_ = new std::vector<pat::Muon>();
  looseElelist_ = new std::vector<pat::Electron>();
  looseEleSIPlist_ = new std::vector<pat::Electron>();
  tightElelist_ = new std::vector<pat::Electron>();
  FSRVec_ = new std::vector<pat::PackedCandidate>();
  looseSIPEleFSRpairVec_ = new std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>();
  looseSIPMuFSRpairVec_ = new std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>();
  tightSIPEleFSRpairVec_ = new std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>();
  tightSIPMuFSRpairVec_ = new std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>();
  tightIsoEleFSRpairVec_ = new std::vector<vhtm::IsoElectron>();
  tightIsoMuFSRpairVec_ = new std::vector<vhtm::IsoMuon>();
  selectedIsoObjectsP4_ = new std::vector<TLorentzVector>();
  looseJetVec_ = new std::vector<pat::Jet>();
}

PhysicsObjectSelector::~PhysicsObjectSelector()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

  delete looseMulist_;
  delete   looseMuSIPlist_;
  delete   tightMulist_; 
  delete   looseElelist_;
  delete   looseEleSIPlist_;
  delete   tightElelist_;
  delete   FSRVec_;
  delete   looseSIPEleFSRpairVec_;
  delete   looseSIPMuFSRpairVec_;
  delete   tightIsoEleFSRpairVec_;
  delete   tightIsoMuFSRpairVec_;
  
}


//
// member functions
//
void
PhysicsObjectSelector::muonSelector(const edm::Handle<pat::MuonCollection>& muons, const reco::Vertex& vit)
{
  looseMulist_->clear();
  looseMuSIPlist_->clear();
  std::vector<pat::PackedCandidate> pftemp;
  for( const auto& mu : *muons ) {
    if(mu.pt() <= 5.) continue;
    if(std::fabs(mu.eta()) > 2.4) continue;
    reco::TrackRef tk = mu.muonBestTrack();
    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    dxyWrtPV = tk->dxy(vit.position());
    dzWrtPV  = tk->dz(vit.position());
    
    if(std::fabs(dxyWrtPV) >= 0.5 )      continue;
    if(std::fabs(dzWrtPV) >= 1.)         continue;
    bool quality = (mu.isGlobalMuon() || ( mu.isTrackerMuon() && mu.numberOfMatches() >= 0)) && mu.muonBestTrackType()!=2 ;
    if(!quality) continue;
    looseMulist_->push_back(mu);
    if(std::fabs(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D)) >= 4.)    continue;
    looseMuSIPlist_->push_back(mu);
    //make a dummy entry here; Filled in fsr block
    looseSIPMuFSRpairVec_->push_back({mu,pftemp});
    bool highPtid = false;
    highPtid = isTrackerHighPt(mu,vit) && mu.pt() > 200;
    bool isTight = mu.isPFMuon() || highPtid;
    if(!isTight)    continue; 
    tightMulist_->push_back(mu);
  }
}

bool PhysicsObjectSelector::isTrackerHighPt(const pat::Muon & mu, const reco::Vertex & primaryVertex) {
  return ( mu.numberOfMatchedStations() > 1 
	   && (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt()) < 0.3 
	   && std::abs(mu.muonBestTrack()->dxy(primaryVertex.position())) < 0.2 
	   && std::abs(mu.muonBestTrack()->dz(primaryVertex.position())) < 0.5 
	   && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 
	   && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 );  
}

//Electron Selection
void PhysicsObjectSelector::electronSelector(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& vit) {
  looseElelist_->clear();
  looseEleSIPlist_->clear();
  tightElelist_->clear();  
  std::vector<pat::PackedCandidate> pftemp;
  
  for (const auto& ele: *electrons) {
    
    bool hasGsfTrack = ele.gsfTrack().isNonnull() ? true : false;
    if (ele.pt() <= 7 ) continue;
    if (fabs(ele.eta()) >= 2.5) continue;
    
    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    if (hasGsfTrack) {
      reco::GsfTrackRef tk = ele.gsfTrack();
      dxyWrtPV = tk->dxy(vit.position());
      dzWrtPV  = tk->dz(vit.position());
    }
    if(dxyWrtPV >= 0.5) continue;
    if(dzWrtPV >= 1.) continue;
    looseElelist_->push_back(ele);
    if(std::fabs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D)) >= 4.)    continue;
    if(!leptonCrosscleaned(ele))   continue;
    looseEleSIPlist_->push_back(ele);
    looseSIPEleFSRpairVec_->push_back({ele,pftemp});
    if(HZZ4lUtil::passBDT(std::fabs(ele.superCluster()->eta()), ele.pt(), ele.userFloat("electronBDT")))
      tightElelist_->push_back(ele);
  }
}

bool PhysicsObjectSelector::leptonCrosscleaned(const pat::Electron& ele) {
  bool flag = true;
  if(verbosity_)  std::cout << "Entering Cross clean" << std::endl;
  for (const auto& mu: *looseMuSIPlist_) {
    if (HZZ4lUtil::getP4(ele).DeltaR(HZZ4lUtil::getP4(mu)) < 0.05) {
      flag = false;
      break;
    }
  }
  return flag;
}

void PhysicsObjectSelector::fsrSelector(const edm::Handle<pat::PackedCandidateCollection>& pfs, const reco::Vertex& vit) {
  //FSRVec_->clear();
  //tightSIPEleFSRpairVec_->clear();
  //tightSIPMuFSRpairVec_->clear();
  if(verbosity_)   std::cout << "Loose Muon/Electron size=" << looseSIPMuFSRpairVec_->size() << "/" << looseSIPEleFSRpairVec_->size() << std::endl;
  for (auto& pfcand: *pfs) {
    if (pfcand.pdgId() != 22 || pfcand.pt() <= 2. || std::fabs(pfcand.eta()) >= 2.4) continue;
    //calculate  Isolation for cone size 30
    std::vector<double> isotemp;   
    LeptonIsoCalculator::calcIsoFromPF(pfcand, pfs, 0.30, isotemp);
    double iso = isotemp.at(0) + isotemp.at(2) + isotemp.at(3) + isotemp.at(4);
    double reliso = iso/pfcand.pt();
    if(reliso >= 1.8)       continue;
    if (!passedSuperClusterVeto(pfcand)) continue;
    if(verbosity_)   std::cout << "Photon pre-selection done" << std::endl;
    //Photons are associated to the closest lepton in the event among all those passing loose ID + SIP cut.
    //Discard photons that do not satisfy the cuts ΔR(γ,l)/ETγ2 < 0.012, and ΔR(γ,l)<0.5 
    int muindx = -1, elindx = -1;
    double dRmin = findClosestLepton(pfcand, muindx, elindx);
    if(verbosity_)   std::cout << "Check Closest lepton done" << std::endl;    
    if (muindx < 0 && elindx < 0) continue;
    
    TLorentzVector pfcandP4 = HZZ4lUtil::getP4(pfcand); 
    double dRovEt2 = dRmin/pfcandP4.Et2();
    if (dRovEt2 >= 0.012 || dRmin >= 0.5)       continue; 
    if(verbosity_)   std::cout << "FSR deltaR done" << std::endl;
    // If more than one photon is associated to the same lepton, the lowest dR(pho,l)/ET_pho^2 is selected.
    
    if (elindx > -1) {  // check electron index first
      if (looseSIPEleFSRpairVec_->at(elindx).second.empty()) 
	looseSIPEleFSRpairVec_->at(elindx).second.push_back(pfcand);
      else {
	TLorentzVector prephoP4 = HZZ4lUtil::getP4(looseSIPEleFSRpairVec_->at(elindx).second.at(0));
	TLorentzVector eleP4 = HZZ4lUtil::getP4(looseSIPEleFSRpairVec_->at(elindx).first);
	if (eleP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
	  looseSIPEleFSRpairVec_->at(elindx).second.clear();
	  looseSIPEleFSRpairVec_->at(elindx).second.push_back(pfcand);
	}
      }
    }
    else if (muindx > -1) {
      if (looseSIPMuFSRpairVec_->at(muindx).second.empty())
	looseSIPMuFSRpairVec_->at(muindx).second.push_back(pfcand);
      else {
	TLorentzVector prephoP4 = HZZ4lUtil::getP4(looseSIPMuFSRpairVec_->at(muindx).second.at(0));
	TLorentzVector muP4 = HZZ4lUtil::getP4(looseSIPMuFSRpairVec_->at(muindx).first);
	if (muP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
	  looseSIPMuFSRpairVec_->at(muindx).second.clear();
	  looseSIPMuFSRpairVec_->at(muindx).second.push_back(pfcand);
	}
	
      }
    }
  }//end loop over pf cand
  for(auto& efsr: *looseSIPEleFSRpairVec_) {
    if(!efsr.second.empty()) FSRVec_->push_back(efsr.second.at(0));
    if(HZZ4lUtil::passBDT(std::fabs(efsr.first.superCluster()->eta()), efsr.first.pt(), efsr.first.userFloat("electronBDT")))
      tightSIPEleFSRpairVec_->push_back(efsr);
  }
  for(auto& mufsr: *looseSIPMuFSRpairVec_) {
    if(!mufsr.second.empty()) FSRVec_->push_back(mufsr.second.at(0));
    bool highPtid = false;
    highPtid = isTrackerHighPt(mufsr.first,vit) && mufsr.first.pt() > 200;
    if( mufsr.first.isPFMuon() || highPtid )
      tightSIPMuFSRpairVec_->push_back(mufsr);
  }
}

bool PhysicsObjectSelector::passedSuperClusterVeto(const pat::PackedCandidate& pfcand, bool verbose) {
  // Supercluster veto by PF reference: veto all the PF candidates used in the PF cluster, 
  // as returned by the method electron.associatedPackedPFCandidates()
  bool passedVeto = true;
  if (verbose && looseEleSIPlist_->size())
    std::cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << std::endl;
  TLorentzVector pfcandP4 = HZZ4lUtil::getP4(pfcand);
  for (const auto& ele: *looseEleSIPlist_) {
    const auto& pfref = ele.associatedPackedPFCandidates();
    for(edm::RefVector<pat::PackedCandidateCollection>::const_iterator it = pfref.begin();  it!=pfref.end(); it++) {
      TLorentzVector p4ref;
      p4ref.SetPtEtaPhiE((*it)->pt(),(*it)->eta(),(*it)->phi(),(*it)->energy());
      if(std::fabs(p4ref.Pt() -pfcandP4.Pt()) < 1e-10 && std::fabs(p4ref.Eta() -pfcandP4.Eta()) < 1e-10 &&  std::fabs(p4ref.Phi() -pfcandP4.Phi()) < 1e-10) {
        passedVeto = false;
        break;
      }
    }
  }
  return passedVeto;
}  

double PhysicsObjectSelector::findClosestLepton(const pat::PackedCandidate& pfPho, int& muindx, int& elindx) {
  TLorentzVector phoP4 = HZZ4lUtil::getP4(pfPho);
  double dRmin = 999.;
  muindx = -1;
  // First consider loose muons
  for (unsigned int i = 0; i < looseSIPMuFSRpairVec_->size(); ++i) {
    const pat::Muon& mu = looseSIPMuFSRpairVec_->at(i).first;
    TLorentzVector muP4 = HZZ4lUtil::getP4(mu);
    double dR = muP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      muindx = i;
    }
  }
  if(verbosity_)   std::cout << "Matching to muons done" << std::endl;
  // Then consider loose electron
  elindx = -1;
  for (unsigned int i = 0; i < looseSIPEleFSRpairVec_->size(); ++i) {
    const pat::Electron& ele = looseSIPEleFSRpairVec_->at(i).first;
    TLorentzVector eleP4 = HZZ4lUtil::getP4(ele);
    double dR = eleP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      elindx = i;
    }
  }
  if(verbosity_)   std::cout << "Matching to electrons done" << std::endl;
  return dRmin;
}

void PhysicsObjectSelector::findIsolatedleptons(const double evrho) {

  for(auto& efsr : *tightSIPEleFSRpairVec_) {
    double eiso = LeptonIsoCalculator::geteleReliso(efsr.first, FSRVec_, evrho);
    if(eiso >= 0.35)   continue;  
    vhtm::IsoElectron ie;
    ie.ele = efsr.first;
    if(efsr.second.empty())   ie.hasfsr = false;
    else {
      ie.hasfsr = true;
      ie.fsr = efsr.second.at(0); 
    }
    ie.relIso = eiso;
    tightIsoEleFSRpairVec_->push_back(ie);
    //fill selobjp4 vector for isolated electrons
    selectedIsoObjectsP4_->push_back(HZZ4lUtil::getP4(ie.ele));
  }
  
  for(auto& mfsr : *tightSIPMuFSRpairVec_) {
    double miso = LeptonIsoCalculator::computeMuonReliso(mfsr.first, FSRVec_);
    if(miso >= 0.35)   continue;  
    vhtm::IsoMuon me;
    me.mu = mfsr.first;
    if(mfsr.second.empty())   me.hasfsr = false;
    else {
      me.hasfsr = true;
      me.fsr = mfsr.second.at(0); 
    }
    me.relIso = miso;
    tightIsoMuFSRpairVec_->push_back(me);
    //fill selobjp4 vector for isolated muons
    selectedIsoObjectsP4_->push_back(HZZ4lUtil::getP4(me.mu));
  }
  
  for(auto& fsr : *FSRVec_) {
    //fill selobjp4 vector for all fsr photons
    selectedIsoObjectsP4_->push_back(HZZ4lUtil::getP4(fsr));
  }
  
  
}

void PhysicsObjectSelector::jetSelector(const edm::Handle<pat::JetCollection>& jets) {
  for(auto& j: *jets) {
    if(j.pt() <= 30)       continue;
    if(std::fabs(j.eta()) >= 4.7)   continue;
    if( !jetLeptonCleaning(HZZ4lUtil::getP4(j)))      continue;
    //if (!isLooseJet(j))    continue;
    looseJetVec_->push_back(j);
  }
}

bool PhysicsObjectSelector::jetLeptonCleaning(const TLorentzVector& jetP4, const double dR) 
{
  for(auto& lp4: *selectedIsoObjectsP4_) {
    if(lp4.DeltaR(jetP4) <= dR)      return false;
  }
  return true;
}

bool PhysicsObjectSelector::isLooseJet(const pat::Jet& jet) {
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

void PhysicsObjectSelector::selectObjects(const edm::Handle<pat::MuonCollection>& muons, 
                                          const edm::Handle<pat::ElectronCollection>& electrons, 
                                          const edm::Handle<pat::PackedCandidateCollection>& pfs, 
                                          const reco::Vertex& vit, 
                                          const double evrho) {
  // order of execution is crucial!
  // muonSelector must precede electronSelector
  muonSelector(muons, vit);
  if(verbosity_)   std::cout << "Muon Selection Done!!" << std::endl;  
  electronSelector(electrons, vit);
  if(verbosity_)   std::cout << "Electron Selection Done!!" << std::endl;  
  // After electrons and muons are found, find photons
  fsrSelector(pfs, vit);
  if(verbosity_)   std::cout << "FSR Selection Done!!" << std::endl;
  // Now find isolated leptons after removing contribution from the associated photons
  findIsolatedleptons(evrho);
  if(verbosity_)   std::cout << "Isolated Leptons Done!!" << std::endl;
  // Jets
  //jetSelector(wt);
}
