#ifndef __AnalysisSpace_OnlineAnalysis_LeptonIsoCalculator_h
#define __AnalysisSpace_OnlineAnalysis_LeptonIsoCalculator_h
// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include "TROOT.h"
// user include files

#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//Root
#include "TLorentzVector.h"

namespace LeptonIsoCalculator {

  double pfiso(const pat::Electron& ele, double eventRho, double fsrPhotonEtSum);
 
  double pfiso(const pat::Muon& mu, double fsrPhotonEtSum);

  double computeMuonReliso(const pat::Muon& mu, const std::vector<pat::PackedCandidate>* fsrVec, 
                           double vetoCone=0.01, double isoCone=0.3, bool verbose=false);
  
  double computeEleReliso(const pat::Electron& ele, const std::vector<pat::PackedCandidate>* fsrVec, const double eventRho, 
                   double vetoCone=0.08, double isoCone=0.3, bool verbose=false);
  double geteleReliso(const pat::Electron& ele, const std::vector<pat::PackedCandidate>* fsrVec, const double eventRho);
  void calcIsoFromPF(const pat::PackedCandidate& v, const edm::Handle<pat::PackedCandidateCollection>& pfs, double cone, std::vector<double>& iso);

}
#endif
