#ifndef __AnalysisSpace_OnlineAnalysis_HZZ4lUtil_h
#define __AnalysisSpace_OnlineAnalysis_HZZ4lUtil_h
// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
//Root
#include "TLorentzVector.h"
#include "TROOT.h"
// user include files

#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


namespace HZZ4lUtil {
  
  const double MZnominal = 91.1876;
  enum ZType {
    mumu = 0, ee, wrong
  };
  enum ZZType {
    mmmm = 0, eeee, eemm, mmee, unknown
  };
  
  template <class T>
    class PtComparatorTL {
  public:
    bool operator()(const T &a, const T &b) const {
      return a.Pt() > b.Pt();
    }
  };
  
  class ZZmasComparator {
  public:
    bool operator()(const vhtm::ZZcandidate &a, const vhtm::ZZcandidate &b) const {
      return a.m4l > b.m4l;
    }
  };

  class ZZkdComparator {
  public:
    bool operator()(const vhtm::ZZcandidate &a, const vhtm::ZZcandidate &b) const {
      return a.Dbkgkin() > b.Dbkgkin();
    }
  };

  class ZZz1massComparator {
  public:
    bool operator()(const vhtm::ZZcandidate &a, const vhtm::ZZcandidate &b) const {
      return std::fabs(a.mZ1 - MZnominal) < std::fabs(b.mZ1 - MZnominal);
    }
  };
  
  template<typename T1, typename T2>
    bool sameFlavourZPair(const T1& Z1Cand, const T2& Z2Cand);
  
  double getEleRhoEffectiveArea(double etax);
  
  bool passBDT(const double fSCeta, const double pt, const double BDT);
  
  //void syncDumper(unsigned long int run, unsigned long int lumi, unsigned long int event, const vhtm::ZZcandidate& ZZ, int nJets,
  //		  double jet1Pt, double jet2Pt, const std::map<std::string,double>& kd, const int category, 
  //		  const double m4lrefit, const double m4lrefiterr, const double weight, std::ostream& os);
  void syncDumper(unsigned long int run, unsigned long int lumi, unsigned long int event, const vhtm::ZZcandidate& ZZ, std::ostream& os);

  void calcIsoFromPF(const pat::PackedCandidate& v, 
                     const edm::Handle<pat::PackedCandidateCollection>& pfs, 
                     double cone, std::vector<double>& iso);
  
  
  template<class T>
    TLorentzVector getP4(const T& pf) {
    TLorentzVector lv;
    lv.SetPtEtaPhiE(pf.pt(), pf.eta(), pf.phi(), pf.energy());
    return lv;
  }
}
#endif
