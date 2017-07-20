#ifndef __AnalysisSpace_OnlineAnalysis_PhysicsObjectSelector_h
#define __AnalysisSpace_OnlineAnalysis_PhysicsObjectSelector_h
// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include "TROOT.h"
// user include files

#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//Root
#include "TLorentzVector.h"

class PhysicsObjectSelector {
  public:
    PhysicsObjectSelector(const bool verbose = false);
    ~PhysicsObjectSelector();
    
    void muonSelector(const edm::Handle<pat::MuonCollection>& muons, const reco::Vertex& vit);
    bool isTrackerHighPt(const pat::Muon & mu, const reco::Vertex & primaryVertex);

    void electronSelector(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& vit);
    bool leptonCrosscleaned(const pat::Electron& ele);
    
    void fsrSelector(const edm::Handle<pat::PackedCandidateCollection>& pfs, const reco::Vertex& vit);
    bool passedSuperClusterVeto(const pat::PackedCandidate& pfcand, bool verbose = false);
    double findClosestLepton(const pat::PackedCandidate& pfPho, int& muindx, int& elindx);
    
    void findIsolatedleptons(const double evrho);

    void jetSelector(const edm::Handle<pat::JetCollection>& jets);
    bool jetLeptonCleaning(const TLorentzVector& jetP4, const double dR = 0.4); 
    bool isLooseJet(const pat::Jet& jet);

    void selectObjects(const edm::Handle<pat::MuonCollection>& muons, 
                       const edm::Handle<pat::ElectronCollection>& electrons, 
                       const edm::Handle<pat::PackedCandidateCollection>& pfs, 
                       const edm::Handle<pat::JetCollection>& jets, 
                       const reco::Vertex& vit, const double evrho);
    
    //access functions
    std::vector<pat::Muon>* looseMulist() { return looseMulist_; }
    std::vector<pat::Muon>* looseMuSIPlist() { return looseMuSIPlist_; }
    std::vector<pat::Muon>* tightMulist() { return tightMulist_; }
    std::vector<pat::Electron>* looseElelist() { return looseElelist_; }
    std::vector<pat::Electron>* looseEleSIPlist() { return looseEleSIPlist_; }
    std::vector<pat::Electron>* tightElelist() { return tightElelist_; }
    std::vector<pat::PackedCandidate>* FSRVec() { return FSRVec_; }
    std::vector<std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>* looseSIPEleFSRpairVec() { return looseSIPEleFSRpairVec_; }
    std::vector<std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>* looseSIPMuFSRpairVec() { return looseSIPMuFSRpairVec_; }
    std::vector<std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>* tightSIPEleFSRpairVec() { return tightSIPEleFSRpairVec_; }
    std::vector<std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>* tightSIPMuFSRpairVec() { return tightSIPMuFSRpairVec_; }
    std::vector<vhtm::IsoElectron>* tightIsoEleFSRpairVec() { return tightIsoEleFSRpairVec_; }
    std::vector<vhtm::IsoMuon>* tightIsoMuFSRpairVec() { return tightIsoMuFSRpairVec_; }
    std::vector<pat::Jet>* looseJetVec() { return looseJetVec_; }
    void printObjects(std::ostream& os = std::cout);
    void printMuonInfo(const pat::Muon&, std::ostream& os =std::cout);
    void printElectronInfo(const pat::Electron& ele, std::ostream& os = std::cout);
  private:
    bool verbosity_;
    std::vector<pat::Muon>* looseMulist_;
    std::vector<pat::Muon>* looseMuSIPlist_;
    std::vector<pat::Muon>* tightMulist_;
    std::vector<pat::Electron>* looseElelist_;
    std::vector<pat::Electron>* looseEleSIPlist_;
    std::vector<pat::Electron>* tightElelist_;
    std::vector<pat::PackedCandidate>* FSRVec_;
    std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>* looseSIPEleFSRpairVec_;
    std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>* looseSIPMuFSRpairVec_;
    std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>* tightSIPEleFSRpairVec_;
    std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>* tightSIPMuFSRpairVec_;
    std::vector<vhtm::IsoElectron>* tightIsoEleFSRpairVec_;
    std::vector<vhtm::IsoMuon>* tightIsoMuFSRpairVec_;
    std::vector<TLorentzVector>* selectedIsoObjectsP4_;//p4 of selected iso objects + fsr
    std::vector<pat::Jet>* looseJetVec_;
    reco::Vertex pVtx_;
    //cut maps
    
};
#endif
