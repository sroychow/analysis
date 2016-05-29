#include <iostream>
#include <algorithm>
#include <utility> 
#include <typeinfo>
#include <iomanip>
#include "TTree.h"
#include "TVector2.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"

namespace HZZ4lUtil {

template<typename T1, typename T2>
 bool sameFlavourZPair( const T1& Z1Cand, const T2& Z2Cand ) {
   if( typeid(Z1Cand) == typeid(Z2Cand) )   return true;
   return false;
}

double getEleRhoEffectiveArea(double etax) {
  double area=0.;
  double eta = std::fabs(etax);
  if (eta >= 0.0 && eta < 0.8) area = 0.1830;
  else if (eta >= 0.8 && eta < 1.3) area = 0.1734;
  else if (eta >= 1.3 && eta < 2.0) area = 0.1077;
  else if (eta >= 2.0 && eta < 2.2) area = 0.1565;
  else if (eta >= 2.2) area = 0.2680;
  return area;
}


double pfiso(const pat::Electron& ele, double eventRho, double fsrPhotonEtSum) {
  return (ele.pfIsolationVariables().sumChargedHadronPt 
          + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt
                       + ele.pfIsolationVariables().sumPhotonEt 
                       - fsrPhotonEtSum - getEleRhoEffectiveArea(std::fabs(ele.eta())) * eventRho));
}

double pfiso(const pat::Muon& mu, double fsrPhotonEtSum) {
  return (mu.pfIsolationR03().sumChargedHadronPt 
          + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt
                       + mu.pfIsolationR03().sumPhotonEt
                       - fsrPhotonEtSum - 0.5 * mu.pfIsolationR03().sumPUPt));
}

  // Isolation in new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone of all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons

  double computeElectronReliso(const pat::Electron& ele, const edm::Handle<pat::PackedCandidateCollection>& fsrVec, double eventRho,
                               double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4 = getP4(ele);
    for (const auto& v: *fsrVec) {
      TLorentzVector fsrP4 = getP4(v);
      double dR = lP4.DeltaR(fsrP4);
      if ((std::fabs(ele.superCluster()->eta()) < 1.479 || dR > vetoCone) && dR < isoCone)
        phoEtSum += fsrP4.Et();
    }
    double iso = pfiso(ele, eventRho, phoEtSum);
    /*
    if (verbose) {
      cout << "electron isolation: " << endl;
      cout << "      iso    lepPt  chHadPt neuHadEt photonEt    fsrEt  effArea eventRho" << endl;
      cout << setprecision(3) 
           << setw(9) << iso 
           << setw(9) << lP4.Pt()
           << setw(9) << ele.chargedHadronIso
           << setw(9) << ele.neutralHadronIso 
           << setw(9) << ele.photonIso
           << setw(9) << phoEtSum
           << setw(9) << getEleRhoEffectiveArea(std::fabs(lP4.Eta()))
           << setw(9) << eventRho
           << endl;
    }*/
    return iso/lP4.Pt();
  }


  // Isolation with new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // dR > 0.01 for muons
  double computeMuonReliso(const pat::Muon& mu, const edm::Handle<pat::PackedCandidateCollection>& fsrVec, 
                           double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4 = getP4(mu);
    for (const auto& v: *fsrVec) {
      TLorentzVector fsrP4 = getP4(v);
      double dR = lP4.DeltaR(fsrP4);
      if (dR > vetoCone && dR < isoCone)
        phoEtSum += fsrP4.Et();
    }
    double iso = pfiso(mu, phoEtSum);
    /*
    if (verbose) {
      cout << "muon isolation: " << endl;
      cout << "      iso    lepPt  chHadPt neuHadEt photonEt    fsrEt     PUPt" << endl;
      cout << setprecision(3) 
           << setw(9) << iso 
           << setw(9) << lP4.Pt()
           << setw(9) << mu.sumChargedHadronPt
           << setw(9) << mu.sumNeutralHadronEt
           << setw(9) << mu.sumPhotonEt
           << setw(9) << phoEtSum
           << setw(9) << mu.sumPUPt
           << endl;
    }
    */
    return iso/lP4.Pt();
  }

bool passBDT(const double fSCeta, const double pt, const double BDT) {
  bool isBDT = (pt<=10 && ((fSCeta<0.8                  && BDT > -0.265) ||
                           (fSCeta>=0.8 && fSCeta<1.479 && BDT > -0.556) ||
                           (fSCeta>=1.479               && BDT > -0.551)))
            || (pt>10  && ((fSCeta<0.8                  && BDT > -0.072) ||
                           (fSCeta>=0.8 && fSCeta<1.479 && BDT > -0.286) ||
                           (fSCeta>=1.479               && BDT > -0.267)));
  return isBDT;
}

/*
double pfiso(const vhtm::PackedPFCandidate& cand) {
  return (cand.isolationMap.at("c30").at(0) 
          + cand.isolationMap.at("c30").at(2) 
          + cand.isolationMap.at("c30").at(3) 
          + cand.isolationMap.at("c30").at(4));
}
*/
void syncDumper(unsigned long int run, unsigned long int lumi, unsigned long int event, const vhtm::ZZcandidate& ZZ, int nJets,
  	        double jet1Pt, double jet2Pt, const std::map<std::string,double>& kd, const int category, 
                const double m4lrefit, const double m4lrefiterr, const double weight, std::ostream& os) {
    //{runX}:{lumiX}:{eventX}:{mass4l:.2fX}:{mZ1:.2fX}:{mZ2:.2fX}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
    //{njets30:dX}: {jet1pt:.2fX}:{jet2pt:.2fX}:{category}
  os << std::fixed << std::setprecision(2);
  os << run << ":"
     << lumi << ":"
     << event << ":"
     << ZZ.m4l << ":"
     << ZZ.mZ1 << ":"
     << ZZ.mZ2 << ":";
   
  os << std::fixed << std::setprecision(3);
  os << kd.find("D_bkg_kin")->second << ":"
     << kd.find("D_bkg")->second  << ":"
     << kd.find("Dgg10_VAMCFM")->second << ":"
     << kd.find("Djet_VAJHU")->second << ":"
     << kd.find("D_g4")->second << ":";
  os << std::fixed << std::setprecision(2);
  os << nJets <<  ":"
     << jet1Pt << ":"
     << jet2Pt << ":"
     << category << ":"
     << m4lrefit << ":"
     << m4lrefiterr << ":"
     << weight
     << std::endl;
  }

}

