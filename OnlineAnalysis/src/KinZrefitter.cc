#include <iostream>
#include <algorithm>
#include <vector>
#include "AnalysisSpace/OnlineAnalysis/src/KinZrefitter.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

#include <TRandom.h>
KinZrefitter::KinZrefitter(const bool mc, const bool verbose) {
  verbosity_ = verbose;
}

KinZrefitter::~KinZrefitter()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}
//
// member functions
//
void KinZrefitter::doZmassrefit(vhtm::ZZcandidate& ZZcand) {
  KinZfitter* kinZfitter = new KinZfitter(!isMC_);
  vector<reco::Candidate *> selectedLeptons;
  TLorentzVector pL11, pL12, pL21, pL22;
  std::map<unsigned int, TLorentzVector> selectedFsrMap;
  TLorentzVector tempP4;
  tempP4.SetPtEtaPhiE(0.,0.,0.,0.);
  if(ZZcand.flavour == HZZ4lUtil::ZZType::eeee) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep2));
  } else if(ZZcand.flavour == HZZ4lUtil::ZZType::mmmm) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep2));
    //if(ZZcand.Z2mm.lep2hasfsr)   selectedFsrMap[3] = HZZ4lUtil::getP4(ZZcand.Z2mm.fsrl2);   
  } else if(ZZcand.flavour == HZZ4lUtil::ZZType::mmee) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1mm.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2ee.lep2));
  } else if(ZZcand.flavour == HZZ4lUtil::ZZType::eemm) {
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z1ee.lep2));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep1));
    selectedLeptons.push_back(dynamic_cast<reco::Candidate* >(&ZZcand.Z2mm.lep2));
  }
  selectedFsrMap[0] = ZZcand.fsrP4Vec[0];
  selectedFsrMap[1] = ZZcand.fsrP4Vec[1];
  selectedFsrMap[2] = ZZcand.fsrP4Vec[2];
  selectedFsrMap[3] = ZZcand.fsrP4Vec[3];
  kinZfitter->Setup(selectedLeptons, selectedFsrMap);
  kinZfitter->KinRefitZ();
  // refit mass4l
  ZZcand.mass4lREFIT(kinZfitter->GetRefitM4l());
  // four 4-vectors after refitting order by Z1_1,Z1_2,Z2_1,Z2_2
  //vector < TLorentzVector > p4 = kinZfitter->GetRefitP4s();//use later 
  // refitted mass4l error
  ZZcand.mass4lErrREFIT(kinZfitter->GetRefitM4lErrFullCov());
  delete kinZfitter; 
}
