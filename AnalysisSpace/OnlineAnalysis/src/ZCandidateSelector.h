#ifndef __AnalysisSpace_OnlineAnalysis_ZCandidateSelector_h
#define __AnalysisSpace_OnlineAnalysis_ZCandidateSelector_h
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
//Root
#include "TLorentzVector.h"

class ZCandidateSelector {
  public:
    ZCandidateSelector(const bool verbose = false);
    ~ZCandidateSelector();
    void ZmumuSelector(const std::vector<vhtm::IsoMuon>& lepPhotonPairVec);
    void ZeeSelector(const std::vector<vhtm::IsoElectron>& lepPhotonPairVec);
    void selectZcandidates(const std::vector<vhtm::IsoMuon>& imuons, const std::vector<vhtm::IsoElectron>& ielectrons);
    std::vector<vhtm::Zee>* getZeeVec() { return ZeeVec_; }
    std::vector<vhtm::Zmumu>* getZmumuVec() { return ZmumuVec_; }  
    //alternate implementation
    void zCandidateSelector(const std::vector<pat::Electron>& tightElelist, const std::vector<pat::Muon>& tightMulist, std::vector<vhtm::Zcandidate>& zcand);
  private:
    bool verbosity_;
    std::vector<vhtm::Zee>* ZeeVec_;
    std::vector<vhtm::Zmumu>* ZmumuVec_;  
};
#endif
