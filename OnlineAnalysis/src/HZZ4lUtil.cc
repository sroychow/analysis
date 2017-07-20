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
  /*older version
    if (eta >= 0.0 && eta < 0.8) area = 0.1830;
    else if (eta >= 0.8 && eta < 1.3) area = 0.1734;
    else if (eta >= 1.3 && eta < 2.0) area = 0.1077;
    else if (eta >= 2.0 && eta < 2.2) area = 0.1565;
    else if (eta >= 2.2) area = 0.2680;
  */
  if( eta >= 0.0000       &&  eta < 1.0000 )  area = 0.1752;
  else if( eta >= 1.0000  &&  eta < 1.4790 )  area = 0.1862;
  else if( eta >= 1.4790  &&  eta < 2.0000 )  area = 0.1411;
  else if( eta >= 2.0000  &&  eta < 2.2000 )  area = 0.1534;
  else if( eta >= 2.2000  &&  eta < 2.3000 )  area = 0.1903;
  else if( eta >= 2.3000  &&  eta < 2.4000 )  area = 0.2243;
  else if( eta >= 2.4000  &&  eta < 5.0000 )  area = 0.2687;
  return area;
}
  
  
  bool passBDT(const double fSCeta, const double pt, const double BDT) {
    /*Pre-Ichep 16 working point
      bool isBDT = (pt<=10 && ((fSCeta<0.8                  && BDT > -0.265) ||
      (fSCeta>=0.8 && fSCeta<1.479 && BDT > -0.556) ||
      (fSCeta>=1.479               && BDT > -0.551)))
      || (pt>10  && ((fSCeta<0.8                  && BDT > -0.072) ||
      (fSCeta>=0.8 && fSCeta<1.479 && BDT > -0.286) ||
      (fSCeta>=1.479               && BDT > -0.267)));
      
    */
    bool isBDT = (pt<=10. && ((fSCeta<0.8                  && BDT > -0.211) ||
			      (fSCeta>=0.8 && fSCeta<1.479 && BDT > -0.396) ||
			      (fSCeta>=1.479               && BDT > -0.215))) 
      || (pt>10.  && ((fSCeta<0.8                   && BDT > -0.870) ||
		      (fSCeta>=0.8 && fSCeta<1.479  && BDT > -0.838) || 
		      (fSCeta>=1.479                && BDT > -0.763)));
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
 // void syncDumper(unsigned long int run, unsigned long int lumi, unsigned long int event, const vhtm::ZZcandidate& ZZ, int nJets,
 //		  double jet1Pt, double jet2Pt, const std::map<std::string,double>& kd, const int category, 
 //		  const double m4lrefit, const double m4lrefiterr, const double weight, std::ostream& os) {
  void syncDumper(unsigned long int run, unsigned long int lumi, unsigned long int event, 
                  const vhtm::ZZcandidate& ZZ, const std::vector<pat::Jet>& cleanJets, std::ostream& os) {
  //{run}:{lumi}:{event}:{mass4l:.2f}:{mZ1:.2f}:{mZ2:.2f}::
  //{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{Dkin_HJJ^VBF:.3f}:{D_0-:.3f}:{Dkin_HJ^VBF-1:.3f}:    
  //{Dkin_HJJ^WH-h:.3f}:{Dkin_HJJ^ZH-h:.3f}:
  //{njets30:d}:{jet1pt:.2f}:{jet2pt:.2f}:{jet1qgl:.3f}:{jet2qgl:.3f}:
  //{Dfull_HJJ^VBF:.3f}:{Dfull_HJ^VBF-1:.3f}:
  //{Dfull_HJJ^WH-h:.3f}:{Dfull_HJJ^ZH-h:.3f}:
  //{category}:
  //{m4lRefit:.2f}:{m4lRefitError:.2f}:{weight:.3f}
    os << std::fixed << std::setprecision(2);
    os << run << ":"
       << lumi << ":"
       << event << ":"
       << ZZ.m4l << ":"
       << ZZ.mZ1 << ":"
       << ZZ.mZ2 << ":";
    os << std::fixed << std::setprecision(3);
    os << ZZ.kdmap.at("D_bkg_kin") << ":"
       << ZZ.kdmap.at("D_bkg") << ":"
       << ZZ.kdmap.at("D_g1g4") << ":"
       << ZZ.kdmap.at("D_VBF") << ":"
       << ZZ.kdmap.at("D_g4") << ":"
       << ZZ.kdmap.at("D_VBF1j") << ":"
       << ZZ.kdmap.at("D_HadWH") << ":"
       << ZZ.kdmap.at("D_HadZH") << ":";
    //printing jet info
    float j1pt = -1., j2pt = -1.;
    float j1qg = -1., j2qg = -1.;
    if(cleanJets.size() == 1) {
      j1pt = cleanJets.at(0).pt();
      j1qg = cleanJets.at(0).userFloat("qgLikelihood");
    }
    else if(cleanJets.size() >= 2) {
      j1pt = cleanJets.at(0).pt();
      j2pt = cleanJets.at(1).pt();
      j1qg = cleanJets.at(0).userFloat("qgLikelihood");
      j2qg = cleanJets.at(1).userFloat("qgLikelihood");
    }
    os << std::fixed << std::setprecision(2);
    os << cleanJets.size() << ":"
       << j1pt << ":"
       << j2pt << ":";
    os << std::fixed << std::setprecision(3);
    os << j1qg << ":"
       << j2qg << ":"
       << ZZ.kdmap.at("D_VBF2j_QG") << ":"
       << ZZ.kdmap.at("D_VBF1j_QG") << ":"
       << ZZ.kdmap.at("D_HadWH_QG") << ":"
       << ZZ.kdmap.at("D_HadZH_QG") << ":";
    os << std::fixed << std::setprecision(2);
    os << -1 << ":"; //this is for category
    os << ZZ.mass4lREFIT() << ":";
    os << ZZ.mass4lErrREFIT() << ":";
    //os << ZZ.flavour;
    os << std::endl;
  }
  void calcIsoFromPF(const pat::PackedCandidate& v, 
		     const edm::Handle<pat::PackedCandidateCollection>& pfs, 
		     double cone, std::vector<double>& iso) {
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
    iso.push_back(chargedHadSum);
    iso.push_back(chargedParticleSum);
    iso.push_back(neutralSum);
    iso.push_back(photonSum);
    iso.push_back(pileupSum);
    //std::cout << "Point cEnd" << std::endl;  
  }
}
