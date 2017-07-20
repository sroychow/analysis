#ifndef __AnalysisSpace_OnlineAnalysis_ZZCandidateSelector_h
#define __AnalysisSpace_OnlineAnalysis_ZZCandidateSelector_h
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
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"
//Root
#include "TLorentzVector.h"

class ZZCandidateSelector {
  public:
    ZZCandidateSelector(const bool verbose = false);
    ~ZZCandidateSelector();
    std::vector<vhtm::ZZcandidate>* getZZVec() { return ZZVec_; }
    void selectZZcandidates(std::vector<vhtm::Zee>&,std::vector<vhtm::Zmumu>&, std::vector<HZZ4lUtil::zzFail>&);
    template <typename T1,typename T2>
    int ZZSelector(const T1& Zcand1, const T2& Zcand2, int& whichZ1cand,double& m4l, bool sameFlavour); 
  private:
    bool verbosity_;
    std::vector<vhtm::ZZcandidate>* ZZVec_;
};
#endif
