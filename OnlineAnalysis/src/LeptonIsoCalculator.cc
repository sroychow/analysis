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

#include "AnalysisSpace/OnlineAnalysis/src/LeptonIsoCalculator.h"

namespace LeptonIsoCalculator {

double pfiso(const pat::Electron& ele, double eventRho, double fsrPhotonEtSum) {
  return (ele.pfIsolationVariables().sumChargedHadronPt 
          + std::max(0., ele.pfIsolationVariables().sumNeutralHadronEt
                       + ele.pfIsolationVariables().sumPhotonEt 
                       - fsrPhotonEtSum - HZZ4lUtil::getEleRhoEffectiveArea(std::fabs(ele.eta())) * eventRho));
}

double pfiso(const pat::Muon& mu, double fsrPhotonEtSum) {
  return (mu.pfIsolationR03().sumChargedHadronPt 
          + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt
                       + mu.pfIsolationR03().sumPhotonEt
                       - fsrPhotonEtSum - 0.5 * mu.pfIsolationR03().sumPUPt));
}

  // Isolation with new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // dR > 0.01 for muons
  double computeMuonReliso(const pat::Muon& mu, const std::vector<pat::PackedCandidate>* fsrVec, 
                           double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4 = HZZ4lUtil::getP4(mu);
    for (const auto& v: *fsrVec) {
      TLorentzVector fsrP4 = HZZ4lUtil::getP4(v);
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

  double geteleReliso(const pat::Electron& ele, const std::vector<pat::PackedCandidate>* fsrVec, const double eventRho) {
    double vetoCone = 0.08;
    double isoCone = 0.3;
    double phoEtSum = 0.;
    TLorentzVector lP4 = HZZ4lUtil::getP4(ele);
    for (const auto& v: *fsrVec) {
      TLorentzVector fsrP4 = HZZ4lUtil::getP4(v);
      double dR = lP4.DeltaR(fsrP4);
      if(dR >= isoCone)   continue;
      if (std::fabs(ele.superCluster()->eta()) < 1.479 || dR > vetoCone)
        phoEtSum += fsrP4.Et();
    }
    double iso = pfiso(ele, eventRho, phoEtSum);
    return iso/lP4.Pt();
  }
  // Isolation in new FSR recovery scheme
  // For each FSR photon that was selected, exclude that photon from the isolation cone of all leptons in the event 
  // passing loose ID + SIP cut if it was in the isolation cone and outside the isolation veto
  // (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons
  double computeEleReliso(const pat::Electron& ele, const std::vector<pat::PackedCandidate>* fsrVec, const double eventRho,
                               double vetoCone, double isoCone, bool verbose)
  {
    double phoEtSum = 0.;
    TLorentzVector lP4 = HZZ4lUtil::getP4(ele);
    for (const auto& v: *fsrVec) {
      TLorentzVector fsrP4 = HZZ4lUtil::getP4(v);
      double dR = lP4.DeltaR(fsrP4);
      if(dR >= isoCone)   continue;
      if (std::fabs(ele.superCluster()->eta()) < 1.479 || dR > vetoCone)
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

void calcIsoFromPF(const pat::PackedCandidate& v, const edm::Handle<pat::PackedCandidateCollection>& pfs, 
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
  //if (verbosity_) std::cout << "isoValues: (" << chargedHadSum << "," 
  //			    << neutralSum << "," << photonSum << "," 
  //			    << pileupSum << ")" 
  //			    << std::endl;
  iso.push_back(chargedHadSum);
  iso.push_back(chargedParticleSum);
  iso.push_back(neutralSum);
  iso.push_back(photonSum);
  iso.push_back(pileupSum);
  //std::cout << "Point cEnd" << std::endl;  
}


}
