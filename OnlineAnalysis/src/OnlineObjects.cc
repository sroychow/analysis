#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

ClassImp(vhtm::IsoElectron)
ClassImp(vhtm::IsoMuon)
ClassImp(vhtm::Zcandidate)
ClassImp(vhtm::Zmumu)
ClassImp(vhtm::Zee)
ClassImp(vhtm::ZZcandidate)
ClassImp(vhtm::ZZcandidate)
ClassImp(vhtm::SelectedEvent)
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

vhtm::Zcandidate::Zcandidate() {
  l1P4.SetPtEtaPhiE(0.,0.,0.,0);
  l2P4.SetPtEtaPhiE(0.,0.,0.,0);
  l1charge = 0;  
  l2charge = 0;
  mass = -999.;
  zmassdiff = -999.;
  //index of the leptons from tight lepton list
  l1idx = -1;
  l2idx = -1;
  flavour = HZZ4lUtil::ZType::wrong;
}
/*
vhtm::Zcandidate::ZZcandidate() {
  //Z1 Z2 will be set according to HZZ selections
  mass = -999.;
  flavour = HZZ4lUtil::ZZType::unknown;
}
*/
vhtm::Zmumu::Zmumu() :
      lep1hasfsr(false),
      lep2hasfsr(false),
      lep1Iso(999.),
      lep2Iso(999.),
      mass(-1.),
      massdiff(9999.)
{
}
void vhtm::Zmumu::dump(std::ostream& os) {
  os << std::setw(8) << mass << std::setw(8) << massdiff << std::endl;
  os << "Lepton 1:";
  HZZ4lUtil::printP4(lep1, os);
  os << "Lep1 has FSR?=" << std::setw(4) << lep1hasfsr << " Relative Isolation=" << lep1Iso << std::endl;
  if(lep1hasfsr)  {
    os << "Attached FSR P4:-";
    HZZ4lUtil::printP4(fsrl1, os);
  }
  os << "Lepton 2:";
  HZZ4lUtil::printP4(lep2, os);
  os << "Lep2 has FSR?=" << std::setw(4) << lep2hasfsr << " Relative Isolation=" << lep2Iso << std::endl;
  if(lep2hasfsr)  {
    os << "Attached FSR P4:-";
    HZZ4lUtil::printP4(fsrl2, os);
  }
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

void vhtm::Zee::dump(std::ostream& os) {
  os << std::setw(8) << mass << std::setw(8) << massdiff << std::endl;
  os << "Lepton 1:";
  HZZ4lUtil::printP4(lep1, os);
  os << "Lep1 has FSR?=" << std::setw(4) << lep1hasfsr << " Relative Isolation=" << lep1Iso << std::endl;
  if(lep1hasfsr)  {
    os << "Attached FSR P4:-";
    HZZ4lUtil::printP4(fsrl1, os);
  }
  os << "Lepton 2:";
  HZZ4lUtil::printP4(lep2, os);
  os << "Lep2 has FSR?=" << std::setw(4) << lep2hasfsr << " Relative Isolation=" << lep2Iso << std::endl;
  if(lep2hasfsr)  {
    os << "Attached FSR P4:-";
    HZZ4lUtil::printP4(fsrl2, os);
  }
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
    if(Z1ee.lep1hasfsr)   
      fsrP4Vec.push_back(HZZ4lUtil::getP4(Z1ee.fsrl1));
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
  //lepton id vector
  lepcharge.push_back(id11);
  lepcharge.push_back(id12);
  lepcharge.push_back(id21);
  lepcharge.push_back(id22);
  //lepton TLorentzVector with fsr
  lepP4VecWithfsr.push_back(lepP4Vec[0] + fsrP4Vec[0]);
  lepP4VecWithfsr.push_back(lepP4Vec[1] + fsrP4Vec[1]);
  lepP4VecWithfsr.push_back(lepP4Vec[2] + fsrP4Vec[2]);
  lepP4VecWithfsr.push_back(lepP4Vec[3] + fsrP4Vec[3]);
}

void vhtm::ZZcandidate::dump(std::ostream& os) {
 os << "Mass Z1=" << std::setw(8) << mZ1;
 os << "Mass Z2=" << std::setw(8) <<  mZ2 << std::endl; 
 if(flavour == HZZ4lUtil::ZZType::eeee) {
   Z1ee.dump(os);
   Z2ee.dump(os);
 } else if(flavour == HZZ4lUtil::ZZType::mmmm) {
   Z1mm.dump(os);
   Z2mm.dump(os);
 } else if(flavour == HZZ4lUtil::ZZType::mmee) {  
   Z1mm.dump(os);
   Z2ee.dump(os);
 } else if(flavour == HZZ4lUtil::ZZType::eemm) { 
   Z1ee.dump(os);
   Z2mm.dump(os);
 }
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
/*
vhtm::SelectedObjects::SelectedObjects() {
    looseMulist_ = new std::vector<pat::Muon>();
    looseMuSIPlist_ = new std::vector<pat::Muon>();
    tightMulist_ = new std::vector<pat::Muon>();
    looseElelist_ = new std::vector<pat::Electron>();
    looseEleSIPlist_ = new std::vector<pat::Electron>();
    tightElelist_ = new std::vector<pat::Electron>();
    FSRVec_ = new std::vector<pat::PackedCandidate>();
    looseSIPEleFSRpairVec_ = new std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>();
    looseSIPMuFSRpairVec_ = new std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>();
    tightSIPEleFSRpairVec_ = new std::vector< std::pair< pat::Electron, std::vector<pat::PackedCandidate> >>();
    tightSIPMuFSRpairVec_ = new std::vector< std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>();
    tightIsoEleFSRpairVec_ = new std::vector<vhtm::IsoElectron>();
    tightIsoMuFSRpairVec_ = new std::vector<vhtm::IsoMuon>();
    selectedIsoObjectsP4_ = new std::vector<TLorentzVector>();//p4 of selected iso objects + fsr
    looseJetVec_ = new std::vector<pat::Jet>();
}
vhtm::SelectedObjects::~SelectedObjects() {
    delete looseMulist_;
    delete looseMuSIPlist_;
    delete tightMulist_;
    delete looseElelist_;
    delete looseEleSIPlist_;
    delete tightElelist_;
    delete FSRVec_;
    delete looseSIPEleFSRpairVec_;
    delete looseSIPMuFSRpairVec_;
    delete tightSIPEleFSRpairVec_;
    delete tightSIPMuFSRpairVec_;
    delete tightIsoEleFSRpairVec_;
    delete tightIsoMuFSRpairVec_;
    delete selectedIsoObjectsP4_;
    delete looseJetVec_;
}
*/
