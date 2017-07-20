#ifndef __AnalysisSpace_OnlineAnalysis_KinZrefitter_h
#define __AnalysisSpace_OnlineAnalysis_KinZrefitter_h
// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include "TROOT.h"
// user include files
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//Root
#include "TLorentzVector.h"
//HEADERS FOR REFIT
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

class KinZrefitter {
  public:
    KinZrefitter(const bool mc, const bool verbose = false);
    ~KinZrefitter();
    void doZmassrefit(vhtm::ZZcandidate& ZZcand);
 private:
    bool verbosity_;
    bool isMC_;
};
#endif
