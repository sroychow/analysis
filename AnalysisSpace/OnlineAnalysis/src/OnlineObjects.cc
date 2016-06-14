#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

ClassImp(vhtm::IsoElectron)
ClassImp(vhtm::IsoMuon)
ClassImp(vhtm::Zmumu)
ClassImp(vhtm::Zee)
ClassImp(vhtm::ZZcandidate)
ClassImp(vhtm::SelectedEvent)
ClassImp(vhtm::ZtnP)
vhtm::IsoElectron::IsoElectron() :
    hasfsr(false),
    relIso(999.)
{
}

vhtm::IsoMuon::IsoMuon() :
    hasfsr(false),
    relIso(999.)
{
}

vhtm::Zmumu::Zmumu() :
      lep1hasfsr(false),
      lep2hasfsr(false),
      lep1Iso(999.),
      lep2Iso(999.),
      mass(-1.),
      massdiff(9999.)
{
}
//reco::RecoCandidate& vhtm::Zmumu::getlep1(){ return  lep1;}
//reco::RecoCandidate& vhtm::Zmumu::getlep2(){ return  lep2;}

vhtm::Zee::Zee() :
      lep1hasfsr(false),
      lep2hasfsr(false),
      lep1Iso(999.),
      lep2Iso(999.),
      mass(-1.),
      massdiff(9999.)
{
}

vhtm::ZZcandidate::ZZcandidate() :
      flavour(4),
      m4l(-1.),
      mZ1(-1.),
      mZ2(-1.),
      Z1mdiff(-1.),
      Z2lepPtavg(-1.)
{
}
void vhtm::ZZcandidate::setP4Vectors() 
{  
  int id11, id12, id21, id22;
  TLorentzVector tempP4;
  tempP4.SetPtEtaPhiE(0.,0.,0.,0.); 
  if(flavour == HZZ4lUtil::ZZType::eeee) {
    id11 = (Z1ee.lep1.charge() < 0) ? 11 : -11;
    id12 = (Z1ee.lep2.charge() < 0) ? 11 : -11;
    id21 = (Z2ee.lep1.charge() < 0) ? 11 : -11;
    id22 = (Z2ee.lep2.charge() < 0) ? 11 : -11;
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.lep2));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.lep2));
    if(Z1ee.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z1ee.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.fsrl2));
    else fsrP4Vec.push_back(tempP4);
    if(Z2ee.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z2ee.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.fsrl2));
    else fsrP4Vec.push_back(tempP4);
  } else if(flavour == HZZ4lUtil::ZZType::mmmm) {
    id11 = (Z1mm.lep1.charge() < 0) ? 13 : -13;
    id12 = (Z1mm.lep2.charge() < 0) ? 13 : -13;
    id21 = (Z2mm.lep1.charge() < 0) ? 13 : -13;
    id22 = (Z2mm.lep2.charge() < 0) ? 13 : -13;
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.lep2));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.lep2));
    if(Z1mm.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z1mm.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.fsrl2));
    else fsrP4Vec.push_back(tempP4);
    if(Z2mm.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z2mm.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.fsrl2));  
    else fsrP4Vec.push_back(tempP4);
  } else if(flavour == HZZ4lUtil::ZZType::mmee) {
    id11 = (Z1mm.lep1.charge() < 0) ? 13 : -13;
    id12 = (Z1mm.lep2.charge() < 0) ? 13 : -13;
    id21 = (Z2ee.lep1.charge() < 0) ? 11 : -11;
    id22 = (Z2ee.lep2.charge() < 0) ? 11 : -11;
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.lep2));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.lep2));
    if(Z1mm.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z1mm.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1mm.fsrl2));
    else fsrP4Vec.push_back(tempP4);
    if(Z2ee.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z2ee.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2ee.fsrl2));   
    else fsrP4Vec.push_back(tempP4);
  } else if(flavour == HZZ4lUtil::ZZType::eemm) {
    id11 = (Z1ee.lep1.charge() < 0) ? 11 : -11;
    id12 = (Z1ee.lep2.charge() < 0) ? 11 : -11;
    id21 = (Z2mm.lep1.charge() < 0) ? 13 : -13;
    id22 = (Z2mm.lep2.charge() < 0) ? 13 : -13;
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.lep2));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.lep1));
    lepP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.lep2));
    if(Z1ee.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z1ee.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.fsrl2));
    else fsrP4Vec.push_back(tempP4);
    if(Z2mm.lep1hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.fsrl1));
    else fsrP4Vec.push_back(tempP4);
    if(Z2mm.lep2hasfsr)   fsrP4Vec.push_back(HZZ4lUtil::getP4(Z2mm.fsrl2)); 
    else fsrP4Vec.push_back(tempP4);
  }
  //lepton charge vector
  lepcharge.push_back(id11);
  lepcharge.push_back(id12);
  lepcharge.push_back(id21);
  lepcharge.push_back(id22);
}

vhtm::SelectedEvent::SelectedEvent() : 
      run(0),
      lumi(0),
      event(0),
      mass4l(0.),
      mZ1(0.),
      mZ2(0.),
      nJets(0),
      jet1pt(0.),
      jet2pt(0.),
      category(0),
      m4lRefit(0.),
      m4lRefitError(0.),
      weight(0.),
      flavour(4)	
{
      kd["Dgg10_VAMCFM"] = -1.; 
      kd["D_bkg_kin"] = -1.; 
      kd["D_bkg"] = -1.; 
      kd["D_g4"] = -1.; 
      kd["Djet_VAJHU"] = -1.;
}

void vhtm::SelectedEvent::reset() 
{
      run = 0;
      lumi = 0;
      event = 0;
      mass4l = 0.;
      mZ1 = 0.;
      mZ2 = 0.;
      nJets = 0.;
      jet1pt = 0.;
      jet2pt = 0.;
      category = 0;
      m4lRefit = 0.;
      m4lRefitError = 0.;
      weight = 0.;
      kd["Dgg10_VAMCFM"] = -1.; 
      kd["D_bkg_kin"] = -1.; 
      kd["D_bkg"] = -1.; 
      kd["D_g4"] = -1.; 
      kd["Djet_VAJHU"] = -1.;
      flavour = 4;
}

vhtm::ZtnP::ZtnP() 
{

}
  

