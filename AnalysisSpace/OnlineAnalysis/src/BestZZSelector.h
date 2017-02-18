#ifndef __AnalysisSpace_OnlineAnalysis_BestZZSelector_h
#define __AnalysisSpace_OnlineAnalysis_BestZZSelector_h
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

class BestZZSelector {
  public:
    BestZZSelector(const bool verbose = false);
    ~BestZZSelector();
  private:
    bool verbosity_;
};
#endif
