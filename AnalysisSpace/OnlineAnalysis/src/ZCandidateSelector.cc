#include <iostream>
#include <algorithm>

#include "AnalysisSpace/OnlineAnalysis/src/ZCandidateSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

ZCandidateSelector::ZCandidateSelector(const bool verbose) {
  verbosity_ = verbose;
  ZeeVec_ = new std::vector<vhtm::Zee>();
  ZmumuVec_ =new std::vector<vhtm::Zmumu>();
}

ZCandidateSelector::~ZCandidateSelector()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  delete ZeeVec_;
  delete ZmumuVec_;
}
//
// member functions
//
void ZCandidateSelector::ZmumuSelector(const std::vector<vhtm::IsoMuon>& lepPhotonPairVec) {
  for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
    const auto& ip = lepPhotonPairVec[i];
    const TLorentzVector& lep1P4 = HZZ4lUtil::getP4(ip.mu);
    TLorentzVector lep1fsrP4;
    lep1fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
    if(!ip.hasfsr) lep1fsrP4 = HZZ4lUtil::getP4(ip.fsr);    
    for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
      const auto& jp = lepPhotonPairVec[j];

      if (ip.mu.charge() + jp.mu.charge() != 0) continue; 

      const TLorentzVector& lep2P4 = HZZ4lUtil::getP4(jp.mu);
      TLorentzVector lep2fsrP4;
      lep2fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
      if(!jp.hasfsr) lep2fsrP4 = HZZ4lUtil::getP4(jp.fsr);

      double zM = (lep1P4 + lep2P4 + lep1fsrP4 + lep2fsrP4).M();
      if(zM <= 12. && zM >= 120.)     continue;
      vhtm::Zmumu Ztemp;
      Ztemp.lep1 = ip.mu;
      Ztemp.lep2 = jp.mu;
      Ztemp.fsrl1 = ip.fsr;
      Ztemp.fsrl2 = jp.fsr;
      Ztemp.lep1hasfsr = ip.hasfsr;
      Ztemp.lep2hasfsr = jp.hasfsr;
      Ztemp.lep1Iso = ip.relIso;
      Ztemp.lep2Iso = jp.relIso;
      Ztemp.mass = zM;
      Ztemp.massdiff = std::fabs(HZZ4lUtil::MZnominal - zM);
      ZmumuVec_->push_back(Ztemp);
    }
}
    
}
void ZCandidateSelector::ZeeSelector(const std::vector<vhtm::IsoElectron>& lepPhotonPairVec) {
  for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
    const auto& ip = lepPhotonPairVec[i];
    const TLorentzVector& lep1P4 = HZZ4lUtil::getP4(ip.ele);
    TLorentzVector lep1fsrP4;
    lep1fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
    if(!ip.hasfsr) lep1fsrP4 = HZZ4lUtil::getP4(ip.fsr);    
    for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
      const auto& jp = lepPhotonPairVec[j];

      if (ip.ele.charge() + jp.ele.charge() != 0) continue; 

      const TLorentzVector& lep2P4 = HZZ4lUtil::getP4(jp.ele);
      TLorentzVector lep2fsrP4;
      lep2fsrP4.SetPtEtaPhiE(0.,0.,0.,0.);
      if(!jp.hasfsr) lep2fsrP4 = HZZ4lUtil::getP4(jp.fsr);

      double zM = (lep1P4 + lep2P4 + lep1fsrP4 + lep2fsrP4).M();
      if(zM <= 12. && zM >= 120.)     continue;
      vhtm::Zee Ztemp;
      Ztemp.lep1 = ip.ele;
      Ztemp.lep2 = jp.ele;
      Ztemp.fsrl1 = ip.fsr;
      Ztemp.fsrl2 = jp.fsr;
      Ztemp.lep1hasfsr = ip.hasfsr;
      Ztemp.lep2hasfsr = jp.hasfsr;
      Ztemp.lep1Iso = ip.relIso;
      Ztemp.lep2Iso = jp.relIso;
      Ztemp.mass = zM;
      Ztemp.massdiff = std::fabs(HZZ4lUtil::MZnominal - zM);
      ZeeVec_->push_back(Ztemp);
    }
  }  
}

void ZCandidateSelector::selectZcandidates(const std::vector<vhtm::IsoMuon>& imuons, const std::vector<vhtm::IsoElectron>& ielectrons) {
  ZmumuSelector(imuons);
  ZeeSelector(ielectrons);
}
//alternate implementation
void zCandidateSelector(const std::vector<pat::Electron>& tightElelist, const std::vector<pat::Muon>& tightMulist,     
                        std::vector<vhtm::Zcandidate>& zcand) {
  //first loop over muons
  for (unsigned int i = 0; i < tightMulist.size(); ++i) {
    const auto& imu = tightMulist[i];
    if(imu.userInt("hzzpassIso") != 1)    continue;
    const TLorentzVector lep1P4 = *imu.userData<TLorentzVector>("dressedLepP4");//fsr is already added if present
    for (unsigned int j = 0; j < tightMulist.size(); ++j) {
      const auto& jmu = tightMulist[j];
      if(jmu.userInt("hzzpassIso") != 1)    continue;
      if(imu.charge() + jmu.charge() != 0)  continue;//OSSF
      const TLorentzVector lep2P4 = *jmu.userData<TLorentzVector>("dressedLepP4");//fsr is already added if present
      double zM = (lep1P4 + lep2P4).M();
      if(zM <= 12. && zM >= 120.)     continue;
      vhtm::Zcandidate ztemp;
      ztemp.l1P4 = HZZ4lUtil::getP4(imu);
      ztemp.l2P4 = HZZ4lUtil::getP4(jmu);
      ztemp.l1wfsrP4 = lep1P4;
      ztemp.l2wfsrP4 = lep2P4;
      ztemp.l1charge = imu.charge();  
      ztemp.l2charge = jmu.charge();
      ztemp.mass = zM;
      ztemp.zmassdiff = std::fabs(HZZ4lUtil::MZnominal - zM);;
      //index of the leptons from tight lepton list
      ztemp.l1idx = i; 
      ztemp.l2idx = j;
      ztemp.flavour = HZZ4lUtil::ZType::mumu;
      zcand.push_back(ztemp);
    }
  }
  //loop over electrons
  for (unsigned int i = 0; i < tightElelist.size(); ++i) {
    const auto& iele = tightElelist[i];
    if(iele.userInt("hzzpassIso") != 1)    continue;
    const TLorentzVector lep1P4 = *iele.userData<TLorentzVector>("dressedLepP4");//fsr is already added if present
    for (unsigned int j = 0; j < tightElelist.size(); ++j) {
      const auto& jele = tightElelist[j];
      if(jele.userInt("hzzpassIso") != 1)    continue;
      if(iele.charge() + jele.charge() != 0)  continue;//OSSF
      const TLorentzVector lep2P4 = *jele.userData<TLorentzVector>("dressedLepP4");//fsr is already added if present
      double zM = (lep1P4 + lep2P4).M();
      if(zM <= 12. && zM >= 120.)     continue;
      vhtm::Zcandidate ztemp;
      ztemp.l1P4 = HZZ4lUtil::getP4(iele);
      ztemp.l2P4 = HZZ4lUtil::getP4(iele);
      ztemp.l1wfsrP4 = lep1P4;
      ztemp.l2wfsrP4 = lep2P4;
      ztemp.l1charge = iele.charge();  
      ztemp.l2charge = jele.charge();
      ztemp.mass = zM;
      ztemp.zmassdiff = std::fabs(HZZ4lUtil::MZnominal - zM);;
      //index of the leptons from tight lepton list
      ztemp.l1idx = i; 
      ztemp.l2idx = j;
      ztemp.flavour = HZZ4lUtil::ZType::ee;
      zcand.push_back(ztemp);
    }
  }

}


