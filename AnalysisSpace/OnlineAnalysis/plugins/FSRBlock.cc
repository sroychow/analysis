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

#include "AnalysisSpace/OnlineAnalysis/plugins/FSRBlock.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"


//
// constructors and destructor
//
FSRBlock::FSRBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfcandTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfCands",edm::InputTag("packedPFCandidates"))), 
  looseElectronTag_(iConfig.getUntrackedParameter<edm::InputTag>("looseElectronSrc",edm::InputTag("looseSIPElectronVector"))),
  looseMuonTag_(iConfig.getUntrackedParameter<edm::InputTag>("looseMuonSrc",edm::InputTag("looseSIPMuonVector"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_)),
  looseElectronToken_(consumes<pat::ElectronCollection>(looseElectronTag_)),
  looseMuonToken_(consumes<pat::MuonCollection>(looseMuonTag_)),
  selectedFSRColl_(iConfig.getUntrackedParameter<std::string>("selectedFSRColl", "selectedFSRVector")),
  looseSIPEleFSRColl_(iConfig.getUntrackedParameter<std::string>("looseSIPEleFSRColl", "looseEleFSRPairList")),
  looseSIPMuFSRColl_(iConfig.getUntrackedParameter<std::string>("looseSIPMuFSRColl", "looseMuFSRPairList"))
{
  produces<std::vector<pat::PackedCandidate>>(selectedFSRColl_);
  produces< std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >> >(looseSIPEleFSRColl_);
  produces< std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >> >(looseSIPMuFSRColl_);
}


FSRBlock::~FSRBlock()
{
}


// ------------ method called to produce the data  ------------
void
FSRBlock::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //std::cout << "Entering FSRBlock::produce" << std::endl;
  FSRVec_->clear();
  looseEleFSRpairVec_->clear();
  looseMuFSRpairVec_->clear();
  
  //using namespace edm;
  
  edm::Handle<pat::PackedCandidateCollection> pfs;
  bool pf_found = iEvent.getByToken(pfToken_, pfs);
  
  edm::Handle<pat::ElectronCollection> leles;
  bool ele_found = iEvent.getByToken(looseElectronToken_, leles);
  //if(ele_found) std::cout<<"Got Loose Electrons"<<std::endl;
  // else std::cout<<"No Loose Electron Found "<<std::endl;
  
  edm::Handle<pat::MuonCollection> lmus;
  bool mu_found=iEvent.getByToken(looseMuonToken_, lmus);
  //if(mu_found) std::cout<<"Got Loose Muons"<<std::endl;
  //else std::cout<<"No Loose Muons Found"<<std::endl;
  
  bool found = pf_found && (mu_found || ele_found);
  //std::cout << "Point 1" << std::endl;  
  //create a dummy collection
  std::vector<pat::PackedCandidate> pfTempv;
  for (unsigned int i = 0; i < leles->size(); ++i)
    looseEleFSRpairVec_->push_back({leles->at(i),pfTempv});
  
  for (unsigned int i = 0; i < lmus->size(); ++i)
    looseMuFSRpairVec_->push_back({lmus->at(i),pfTempv});
 
  //std::cout << "Point 2" << std::endl;  
  
  /*
    Preselection: pT > 2 GeV, |η| < 2.4, photon PF relative isolation less than 1.8.
    The PF isolation is computed using a cone of 0.3, a threshold of
    0.2 GeV on charged hadrons with a veto cone of 0.0001, and 
    0.5 GeV on neutral hadrons and photons with a veto cone of 0.01, 
    including also the contribution from PU vertices (same radius and threshold as per charged isolation) .
    Supercluster veto by PF reference: veto all the PF candidates used in the PF cluster, 
    as returned by the method electron.associatedPackedPFCandidates() 
   (this is different from Moriond 2016 analysis - for details see Giovani's presentation),
  */
  if (found && pfs.isValid()) {
  //std::cout << "Point 2a" << std::endl;  
    if (leles->size() || lmus->size()) { 
  //std::cout << "Point 2b" << std::endl;  
     for (auto& pfcand: *pfs) {
	
	if (pfcand.pdgId() != 22 || pfcand.pt() <= 2. || std::fabs(pfcand.eta()) >= 2.4) continue;
  //std::cout << "Point 2c" << std::endl;  

          //calculate  Isolation for cone size 30
	std::vector<double> isotemp;   
	calcIsoFromPF(pfcand, pfs, 0.30, isotemp);
	double iso = isotemp.at(0)
	  + isotemp.at(2)
	  + isotemp.at(3)
	  + isotemp.at(4);
	double reliso = iso/pfcand.pt();
	
	if(reliso >= 1.8)       continue;
     //std::cout << "Point 2d" << std::endl;  
	
	if (!passedSuperClusterVeto(pfcand, leles)) continue;
	
	//Photons are associated to the closest lepton in the event among all those passing loose ID + SIP cut.
	//Discard photons that do not satisfy the cuts ΔR(γ,l)/ETγ2 < 0.012, and ΔR(γ,l)<0.5 
	int muindx = -1, elindx = -1;
	double dRmin = findClosestLepton(pfcand, leles, lmus,muindx, elindx);
	
	if (muindx < 0 && elindx < 0) continue;
	
	TLorentzVector pfcandP4 = getP4(pfcand); 
	double dRovEt2 = dRmin/pfcandP4.Et2();
	if (dRovEt2 >= 0.012 || dRmin >= 0.5)       continue; 
	
	// If more than one photon is associated to the same lepton, the lowest dR(pho,l)/ET_pho^2 is selected.

	if (elindx > -1) {  // check electron index first
	  if (looseEleFSRpairVec_->at(elindx).second.empty()) 
	    looseEleFSRpairVec_->at(elindx).second.push_back(pfcand);
	  else {
	    TLorentzVector prephoP4 = getP4(looseEleFSRpairVec_->at(elindx).second.at(0));
	    TLorentzVector eleP4 = getP4(looseEleFSRpairVec_->at(elindx).first);
	    if (eleP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
	      looseEleFSRpairVec_->at(elindx).second.clear();
	      looseEleFSRpairVec_->at(elindx).second.push_back(pfcand);
	    }
	  }
	}
	else if (muindx > -1) {
	  if (looseMuFSRpairVec_->at(muindx).second.empty())
	    looseMuFSRpairVec_->at(muindx).second.push_back(pfcand);
	  else {
	    TLorentzVector prephoP4 = getP4(looseMuFSRpairVec_->at(muindx).second.at(0));
	    TLorentzVector muP4 = getP4(looseMuFSRpairVec_->at(muindx).first);
	    if (muP4.DeltaR(prephoP4)/prephoP4.Et2() > dRovEt2) {
	      looseMuFSRpairVec_->at(muindx).second.clear();
	      looseMuFSRpairVec_->at(muindx).second.push_back(pfcand);
	    }
  //std::cout << "Point 2e" << std::endl;  

	  }
	}
      }
    }
  }//end loop over pf cand
  //std::cout << "Point 3" << std::endl;
  
  for(auto& efsr: *looseEleFSRpairVec_) {
    if(!efsr.second.empty())
      FSRVec_->push_back(efsr.second.at(0));
  }
  for(auto& mfsr: *looseMuFSRpairVec_) {
    if(!mfsr.second.empty())
      FSRVec_->push_back(mfsr.second.at(0));
  }
  //std::cout << "Point 4" << std::endl;

  std::auto_ptr<std::vector<pat::PackedCandidate>> pv1(new std::vector<pat::PackedCandidate>(*FSRVec_));
  iEvent.put(pv1,selectedFSRColl_);
  std::auto_ptr<std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>> pv2(new std::vector< std::pair< pat::Electron, 
                                                                                                                std::vector<pat::PackedCandidate> >>(*looseEleFSRpairVec_));
  iEvent.put(pv2,looseSIPEleFSRColl_);
  std::auto_ptr<std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>> pv3(new std::vector< std::pair< pat::Muon, 
                                                                                                 std::vector<pat::PackedCandidate> >>(*looseMuFSRpairVec_));
  iEvent.put(pv3,looseSIPMuFSRColl_);
  ////std::cout << "Leaving FSRBlock::produce" << std::endl;
}

void FSRBlock::calcIsoFromPF(const pat::PackedCandidate& v, 
			     edm::Handle<pat::PackedCandidateCollection>& pfs, 
			     double cone, std::vector<double>& iso)
{
  // initialize sums
  double chargedHadSum = 0., 
    chargedParticleSum = 0., 
    neutralSum = 0., 
    photonSum = 0., 
    pileupSum  = 0;
  //std::cout << "Point c1" << std::endl;  
  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0; i < v.numberOfSourceCandidatePtrs(); ++i)
    footprint.push_back(v.sourceCandidatePtr(i));
  //std::cout << "Point c2" << std::endl;  
  
  // now loop on pf candidates
  for (unsigned int i = 0; i < pfs->size(); ++i) {
    const pat::PackedCandidate& pf = (*pfs)[i];
    double dRcone = deltaR(v, pf);
    if (dRcone < cone) {
      // pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs, i)) != footprint.end()) continue;
      
      double pt = pf.pt(); 
      int pdg = std::abs(pf.pdgId());
      if (pf.charge() == 0) {
        if (pt > 0.5 && dRcone > 0.01) {
          if (pdg == 22)
            photonSum += pt;
          else 
            neutralSum += pt;
        }
      } 
      else {
        if (pt > 0.2 && dRcone > 0.0001) { 
          if ( pf.vertexRef().isNonnull() && pf.fromPV() >= 2) {
          //if (1) {
	    chargedParticleSum += pt;
            if (pdg != 13 && pdg != 11) chargedHadSum += pt;
          } 
	  else 
            pileupSum += pt;
        }
      }
    }
  //std::cout << "Point c3" << std::endl;  

  }
  if (verbosity_) std::cout << "isoValues: (" << chargedHadSum << "," 
			    << neutralSum << "," << photonSum << "," 
			    << pileupSum << ")" 
			    << std::endl;
  iso.push_back(chargedHadSum);
  iso.push_back(chargedParticleSum);
  iso.push_back(neutralSum);
  iso.push_back(photonSum);
  iso.push_back(pileupSum);
  //std::cout << "Point cEnd" << std::endl;  
}


bool FSRBlock::passedSuperClusterVeto(const pat::PackedCandidate& pfcand, edm::Handle<pat::ElectronCollection>& le, bool verbose) {
  // Supercluster veto by PF reference: veto all the PF candidates used in the PF cluster, 
  // as returned by the method electron.associatedPackedPFCandidates()
  bool passedVeto = true;
  if (verbose && le->size())
    std::cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << std::endl;
  for (const auto& ele: *le) {
    const auto& pfref = ele.associatedPackedPFCandidates();
    for(edm::RefVector<pat::PackedCandidateCollection>::const_iterator it = pfref.begin();  it!=pfref.end(); it++) {
      TLorentzVector p4ref;
      p4ref.SetPtEtaPhiE((*it)->pt(),(*it)->eta(),(*it)->phi(),(*it)->energy());
      if(p4ref == HZZ4lUtil::getP4(pfcand)) {
        passedVeto = false;
        break;
      }
    }
  }
  return passedVeto;
}

double FSRBlock::findClosestLepton(const pat::PackedCandidate& pfPho, edm::Handle<pat::ElectronCollection>& le, 
                                   edm::Handle<pat::MuonCollection>& lmu,int& muindx, int& elindx) {
  TLorentzVector phoP4 = getP4(pfPho);
  double dRmin = 999.;
  muindx = -1;
  // First consider loose muons
  for (unsigned int i = 0; i < lmu->size(); ++i) {
    const pat::Muon& mu = lmu->at(i);
    TLorentzVector muP4 = getP4(mu);
    double dR = muP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      muindx = i;
    }
  }
  // Then consider loose electron
  elindx = -1;
  for (unsigned int i = 0; i < le->size(); ++i) {
    const pat::Electron& ele = le->at(i);
    TLorentzVector eleP4 = getP4(ele);
    double dR = eleP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      elindx = i;
    }
  }
  return dRmin;
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
FSRBlock::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
FSRBlock::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  FSRBlock::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  FSRBlock::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
FSRBlock::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  FSRVec_ = new std::vector<pat::PackedCandidate>();
  looseEleFSRpairVec_= new std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>();
  looseMuFSRpairVec_ = new std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>();
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
FSRBlock::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FSRBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FSRBlock);
