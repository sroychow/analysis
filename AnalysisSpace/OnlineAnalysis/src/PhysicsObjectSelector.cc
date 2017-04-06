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
  //for(auto& mu : *muons ) {
  for(unsigned int i = 0; i<muons->size(); i++) {
    pat::Muon mu = muons->at(i);
    reco::TrackRef tk = mu.muonBestTrack();
    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    dxyWrtPV = tk->dxy(vit.position());
    dzWrtPV  = tk->dz(vit.position());
    //add user floats now to avoid recalculating during object dump
    mu.addUserFloat("hzzdxyWrtPV", dxyWrtPV); 
    mu.addUserFloat("hzzdzWrtPV", dzWrtPV); 
    mu.addUserFloat("hzzSIP", mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D));
    bool highPtid = false;
    highPtid = isTrackerHighPt(mu,vit) && mu.pt() > 200.;
    bool isTight = mu.isPFMuon() || highPtid;
    mu.addUserInt("hzzhighPtid", highPtid);
    mu.addUserInt("hzzisTight", isTight);
    //loose selection
    if(mu.pt() <= 5.) continue;//pt > 5 GeV
    if(std::abs(mu.eta()) >= 2.4) continue;//|eta| < 2.4
    if(std::abs(dxyWrtPV) >= 0.5 )      continue;//dxy < 0.5
    if(std::abs(dzWrtPV) >= 1.)         continue;//dz < 1.
    bool quality = (mu.isGlobalMuon() || ( mu.isTrackerMuon() && mu.numberOfMatches() > 0)) && mu.muonBestTrackType()!=2 ;
    if(!quality) continue;
    looseMulist_->push_back(mu);
    //sip cut
    if(std::abs(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D)) >= 4.)    continue;
    looseMuSIPlist_->push_back(mu);
    looseSIPMuFSRpairVec_->push_back({mu,pftemp});
    //if(!isTight)    continue; 
    //tightMulist_->push_back(mu);
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

////////////Electron Selection//////////////////////////
void PhysicsObjectSelector::electronSelector(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& vit) {
  looseElelist_->clear();
  looseEleSIPlist_->clear();
  tightElelist_->clear();  
  std::vector<pat::PackedCandidate> pftemp;
  
  //for (auto& ele: *electrons) {
  for(unsigned int i = 0; i<electrons->size(); i++) {
    pat::Electron ele = electrons->at(i);
    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    if (ele.gsfTrack().isNonnull()) {
      reco::GsfTrackRef tk = ele.gsfTrack();
      dxyWrtPV = tk->dxy(vit.position());
      dzWrtPV  = tk->dz(vit.position());
    }
    //add user floats now to avoid recalculating during object dump
    ele.addUserFloat("hzzdxyWrtPV", dxyWrtPV); 
    ele.addUserFloat("hzzdzWrtPV", dzWrtPV); 
    ele.addUserFloat("hzzSIP", ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D));
    bool isTight = HZZ4lUtil::passBDT(std::abs(ele.superCluster()->eta()), ele.pt(), ele.userFloat("electronBDT"));
    ele.addUserInt("hzzisTight", isTight);
    //loose selection
    if (ele.pt() <= 7 ) continue;//pt < 7
    if (abs(ele.eta()) >= 2.5) continue;//|eta| < 2.5
    if(std::abs(dxyWrtPV) >= 0.5) continue; //dxy < 5. 
    if(std::abs(dzWrtPV) >= 1.) continue;//dz < 1.
    looseElelist_->push_back(ele);
    //sip cut
    if(std::abs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D)) >= 4.)    continue;
    if(!leptonCrosscleaned(ele))   continue;
    looseEleSIPlist_->push_back(ele);
    looseSIPEleFSRpairVec_->push_back({ele,pftemp});
    //tight cut
    //if(isTight)
    //tightElelist_->push_back(ele);
  }
}
//Remove electrons which are within ΔR(eta,phi)<0.05 of a muon passing tight ID && SIP<4 
bool PhysicsObjectSelector::leptonCrosscleaned(const pat::Electron& ele) {
  bool flag = true;
  if(verbosity_)  std::cout << "Entering Cross clean" << std::endl;
  for (const auto& mu: *tightMulist_) {
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
    if (pfcand.pdgId() != 22 || pfcand.pt() <= 2. || std::abs(pfcand.eta()) >= 2.4) continue;
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
        //lowest deltaR/et2 
        //if fsr exits, replace if dRovEt2 is greater compared to current candidate
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
        //lowest deltaR/et2 
        //if fsr exits, replace if dRovEt2 is greater compared to current candidate
	if (muP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
	  looseSIPMuFSRpairVec_->at(muindx).second.clear();
	  looseSIPMuFSRpairVec_->at(muindx).second.push_back(pfcand);
	}
	
      }
    }
  }//end loop over pf cand
  //fill tight leptons with fsr vector and the fsr vector
  //alternative implementation
  for(unsigned int i = 0; i<looseSIPMuFSRpairVec_->size(); i++){
    auto&  mfsr = looseSIPMuFSRpairVec_->at(i);
    bool hasFsr = false;
    TLorentzVector pfP;
    if(!mfsr.second.empty()) {
      pfP = HZZ4lUtil::getP4(mfsr.second.at(0));
      hasFsr = true;
    }  
    else pfP.SetPtEtaPhiE(0.,0.,0.,0.);
    looseMuSIPlist_->at(i).addUserInt("hasFSR", hasFsr);
    looseMuSIPlist_->at(i).addUserData("fsrP4", pfP);
    pfP += HZZ4lUtil::getP4(mfsr.first);//add pf P4 to lep P4.
    looseMuSIPlist_->at(i).addUserData("dressedLepP4", pfP); 
    if(looseMuSIPlist_->at(i).userInt("hzzisTight")) 
      tightMulist_->push_back(looseMuSIPlist_->at(i));  
  }
  for(unsigned int i = 0; i<looseSIPEleFSRpairVec_->size(); i++){
    auto&  efsr = looseSIPEleFSRpairVec_->at(i);
    bool hasFsr = false;
    TLorentzVector pfP;
    if(!efsr.second.empty()) {
      pfP = HZZ4lUtil::getP4(efsr.second.at(0));
      hasFsr = true;
    }  
    else pfP.SetPtEtaPhiE(0.,0.,0.,0.);
    looseEleSIPlist_->at(i).addUserInt("hasFSR", hasFsr);
    looseEleSIPlist_->at(i).addUserData("fsrP4", pfP);
    pfP += HZZ4lUtil::getP4(efsr.first);//add pf P4 to lep P4.
    looseEleSIPlist_->at(i).addUserData("dressedLepP4", pfP); 
    if(looseEleSIPlist_->at(i).userInt("hzzisTight")) 
      tightElelist_->push_back(looseEleSIPlist_->at(i));  
  }

  //****
  for(auto& efsr: *looseSIPEleFSRpairVec_) {
    if(!efsr.second.empty()) {
      FSRVec_->push_back(efsr.second.at(0));
    }
    if(HZZ4lUtil::passBDT(std::abs(efsr.first.superCluster()->eta()), efsr.first.pt(), efsr.first.userFloat("electronBDT")))
      tightSIPEleFSRpairVec_->push_back(efsr);
  }
  for(auto& mufsr: *looseSIPMuFSRpairVec_) {
    if(!mufsr.second.empty()) FSRVec_->push_back(mufsr.second.at(0));
    bool highPtid = false;
    highPtid = isTrackerHighPt(mufsr.first,vit) && mufsr.first.pt() > 200.;
    if( mufsr.first.isPFMuon() || highPtid )
      tightSIPMuFSRpairVec_->push_back(mufsr);
  }
}

bool PhysicsObjectSelector::passedSuperClusterVeto(const pat::PackedCandidate& pfcand, bool verbose) {
  // Supercluster veto by PF reference: veto all the PF candidates used in the PF cluster, 
  // as returned by the method electron.associatedPackedPFCandidates()
  //looseElelist_
  bool passedVeto = true;
  if (verbose && looseEleSIPlist_->size())
    std::cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << std::endl;
  TLorentzVector pfcandP4 = HZZ4lUtil::getP4(pfcand);
  //for (const auto& ele: *looseEleSIPlist_) {
  for (const auto& ele: *looseElelist_) {
    const auto& pfref = ele.associatedPackedPFCandidates();
    for(edm::RefVector<pat::PackedCandidateCollection>::const_iterator it = pfref.begin();  it!=pfref.end(); it++) {
      TLorentzVector p4ref;
      p4ref.SetPtEtaPhiE((*it)->pt(),(*it)->eta(),(*it)->phi(),(*it)->energy());
      if(std::abs(p4ref.Pt() -pfcandP4.Pt()) < 1e-10 && std::abs(p4ref.Eta() -pfcandP4.Eta()) < 1e-10 &&  std::abs(p4ref.Phi() -pfcandP4.Phi()) < 1e-10) {
        passedVeto = false;
        break;
      }
    }
  }
  return passedVeto;
}  

double PhysicsObjectSelector::findClosestLepton(const pat::PackedCandidate& pfPho, int& muindx, int& elindx) {
  TLorentzVector phoP4 = HZZ4lUtil::getP4(pfPho);
  double dRmin = 0.5;
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
  //alternate implementation
  for(auto& te: *tightElelist_) {
    double eiso = LeptonIsoCalculator::geteleReliso(te, FSRVec_, evrho);
    te.addUserFloat("hzzrelIso", eiso);
    if(eiso < 0.35)  te.addUserInt("hzzpassIso", 1);
    else  te.addUserInt("hzzpassIso", 0);
  }
  for(auto& tm: *tightMulist_) {
    double miso = LeptonIsoCalculator::computeMuonReliso(tm, FSRVec_);
    tm.addUserFloat("hzzrelIso", miso);
    if(miso < 0.35)  tm.addUserInt("hzzpassIso", 1);
    else  tm.addUserInt("hzzpassIso", 0);
  }

  //
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
    if(std::abs(j.eta()) >= 4.7)   continue;
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
  bool centralCut = (std::abs(jet.eta()) <= 2.4) 
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
                                          const edm::Handle<pat::JetCollection>& jets,  
                                          const reco::Vertex& vit, 
                                          const double evrho) {
  pVtx_ = vit;
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
  jetSelector(jets);
  std::sort( (*looseJetVec_).begin(), (*looseJetVec_).end(), HZZ4lUtil::PtComparatorPAT<pat::Jet>() );
}

void PhysicsObjectSelector::printObjects(std::ostream& os) {
  std::cout << "<<<<Vertex Info>>>>" << std::endl;
  os << std::setw(6) << "ndf"
     << std::setw(6) << "z"
     << std::setw(6) << "chi2"
     << std::setw(6) << "rho" 
  << std::endl;
  os << std::setprecision(2)
     << std::setw(8) << pVtx_.ndof()
     << std::setw(8) << pVtx_.position().Z() 
     << std::setw(8) << pVtx_.chi2()
     << std::setw(8) << pVtx_.position().Rho()
  << std::endl;
  std::cout << "<<<<MUON INFO>>>>" << std::endl;
  os  << std::setw(6) << "Charge"
      << std::setw(8) << "pt"
      << std::setw(8) << "eta"
      << std::setw(8) << "phi"
      << std::setw(8) << "energy"
      << std::setw(9) << "dxy"
      << std::setw(9) << "dz"
      << std::setw(9) << "SIP"
      << std::setprecision(1)
      << std::setw(9) << "isGlobal"
      << std::setw(9) << "isTracker"
      << std::setw(9) << "nMatches"
      << std::setw(6) << "isPF"
      << std::setw(6) << "HptId"
      << std::setw(6) << "isTight"
      << std::endl;  
  std::cout << "Print Loose Muons::" << std::endl;
  for(const auto& mu : *looseMulist_)  printMuonInfo(mu, os);
  std::cout << "Print Tight Muons::" << std::endl;
  for(const auto& mu : *tightMulist_)  printMuonInfo(mu, os);
  std::cout << "<<<<ELECTRON INFO>>>>" << std::endl;
  os  << std::setw(6) << "Charge"
      << std::setw(8) << "pt"
      << std::setw(8) << "eta"
      << std::setw(8) << "phi"
      << std::setw(8) << "energy"
      << std::setw(9) << "dxy"
      << std::setw(9) << "dz"
      << std::setw(9) << "SIP"
      << std::setprecision(1)
      << std::setw(8) << "scEta"
      << std::setw(8) << "BDT"
      << std::setw(8) << "passBDT"
      << std::endl;  
  std::cout << "Print Loose Elelist_::" << std::endl;
  for(const auto& ele : *looseElelist_)  printElectronInfo(ele, os);
  std::cout << "Print Tight Electrons::" << std::endl;
  for(const auto& ele : *tightElelist_)  printElectronInfo(ele, os);
  std::cout << "<<<<FSR INFO>>>>" << std::endl;
  std::cout << "Print Selected FSR photons::" << std::endl;  
  for(const auto& fsr : *FSRVec_)   HZZ4lUtil::printP4(fsr, os);
  
  std::cout << "<<<<<<Print Isolated Tight Muons>>>>>>" << std::endl;
  for(const auto& imu : *tightIsoMuFSRpairVec_) {
    printMuonInfo(imu.mu);
    os << std::setprecision(5);
    os << "Has FSR?=" << std::setw(4) << imu.hasfsr << " RelIso" << imu.relIso << std::endl;
    if(imu.hasfsr)  {
      os << "Attached FSR P4:-";
      HZZ4lUtil::printP4(imu.fsr, os);
    }
  }
  std::cout << "<<<<<<Print Isolated Tight Electrons>>>>>>" << std::endl;  for(const auto& iele : *tightIsoEleFSRpairVec_) {
    printElectronInfo(iele.ele);
    os << std::setprecision(5);
    os << "Has FSR?=" << std::setw(4) << iele.hasfsr << " RelIso" << iele.relIso << std::endl;
    if(iele.hasfsr)  {
      os << "Attached FSR P4:-";
      HZZ4lUtil::printP4(iele.fsr, os);
    }
  }
}

//print muon info
void PhysicsObjectSelector::printMuonInfo(const pat::Muon& muon, std::ostream& os) {
    os << std::setprecision(3)
      << std::setw(6) << muon.charge()
      << std::setw(8) << muon.pt()
      << std::setw(8) << muon.eta()
      << std::setw(8) << muon.phi()
      << std::setw(8) << muon.energy()
      << std::setw(9) << std::abs(muon.muonBestTrack()->dxy(pVtx_.position()))
      << std::setw(8) << std::abs(muon.muonBestTrack()->dz(pVtx_.position()))
      << std::setw(8) << std::abs(muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D))
      << std::setprecision(1)
      << std::setw(8) << muon.isGlobalMuon()
      << std::setw(8) << muon.isTrackerMuon()
      << std::setw(8) << muon.numberOfMatches()
      << std::setw(8) << muon.isPFMuon();
      bool highPtid = false;
      highPtid = isTrackerHighPt(muon,pVtx_) && muon.pt() > 200.;
      bool isTight = muon.isPFMuon() || highPtid;
   os << std::setw(8) << highPtid
      << std::setw(8) << isTight
      << std::endl;
}

//print electron info
void PhysicsObjectSelector::printElectronInfo(const pat::Electron& ele, std::ostream& os) {
    os << std::setprecision(3)
      << std::setw(6) << ele.charge()
      << std::setw(8) << ele.pt()
      << std::setw(8) << ele.eta()
      << std::setw(8) << ele.phi()
      << std::setw(8) << ele.energy()
      << std::setw(8) << std::abs(ele.gsfTrack()->dxy(pVtx_.position()))
      << std::setw(8) << std::abs(ele.gsfTrack()->dz(pVtx_.position()))
      << std::setw(8) << std::abs(ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D))
      << std::setw(8) << std::abs( ele.superCluster()->eta() )
      << std::setw(8) << ele.userFloat("electronBDT")
      << std::setprecision(1)
      << std::setw(8) << HZZ4lUtil::passBDT(std::abs(ele.superCluster()->eta()), ele.pt(), ele.userFloat("electronBDT"))
      << std::endl;
}
