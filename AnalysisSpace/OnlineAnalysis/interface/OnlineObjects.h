#ifndef __AnalysisSpace_OnlineAnalysis_OnlineObjects_h
#define __AnalysisSpace_OnlineAnalysis_OnlineObjects_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "TLorentzVector.h"
#include "TObject.h"

namespace vhtm {

class IsoElectron {
  public:
    IsoElectron();
    virtual ~IsoElectron(){}
    pat::Electron ele;
    bool hasfsr;
    pat::PackedCandidate fsr;
    double relIso;
    ClassDef(IsoElectron,1)
};

class IsoMuon {
  public:
    IsoMuon();
    virtual ~IsoMuon(){}
    pat::Muon mu;
    bool hasfsr;
    pat::PackedCandidate fsr;
    double relIso;
    ClassDef(IsoMuon,1)
};

class Zmumu {
  public:
    Zmumu();
    virtual ~Zmumu(){}
    pat::Muon lep1;
    pat::Muon lep2;
    pat::PackedCandidate fsrl1;
    pat::PackedCandidate fsrl2;
    bool lep1hasfsr;
    bool lep2hasfsr;
    double lep1Iso;
    double lep2Iso;
    double mass;
    double massdiff;
    //reco::RecoCandidate& getlep1();
    //reco::RecoCandidate& getlep2();
    ClassDef(Zmumu,1)
};
class Zee {
  public:
    Zee();
    virtual ~Zee(){}
    pat::Electron lep1;
    pat::Electron lep2;
    pat::PackedCandidate fsrl1;
    pat::PackedCandidate fsrl2;
    bool lep1hasfsr;
    bool lep2hasfsr;
    double lep1Iso;
    double lep2Iso;
    double mass;
    double massdiff;
    //reco::RecoCandidate& getlep1();
    //reco::RecoCandidate& getlep2();
    ClassDef(Zee,1)
};
class ZZcandidate {
  public:
    ZZcandidate();
    virtual ~ZZcandidate(){}
    vhtm::Zmumu Z1mm;
    vhtm::Zmumu Z2mm;
    vhtm::Zee Z1ee;
    vhtm::Zee Z2ee;
    int flavour;//0=4mu,1=4e,2=2e2mu,3=2mu2e,4=wrong
    double m4l;
    double mZ1;
    double mZ2;
    double Z1mdiff;
    double Z2lepPtavg;
    std::vector<int> lepcharge;
    std::vector<TLorentzVector> lepP4Vec;
    std::vector<TLorentzVector> fsrP4Vec;
    void setP4Vectors();
    ClassDef(ZZcandidate,1)
};
class SelectedEvent : public TObject {
  public:
    SelectedEvent();
    virtual ~SelectedEvent(){}
    unsigned long int run;
    unsigned long int lumi;
    unsigned long int event;
    double mass4l;
    double mZ1;
    double mZ2;
    std::map<std::string,double>  kd;//D_bkg^kin,D_bkg,D_gg,D_HJJ^VBF,D_0
    int nJets;
    double jet1pt;
    double jet2pt;
    int category;
    double m4lRefit;
    double m4lRefitError;
    double weight;
    int flavour;//0=4mu,1=4e,2=2e2mu,3=2mu2e,4=wrong
    void reset();
    ClassDef(SelectedEvent,1)
};
class ZtnP : public TObject {
   public :
  ZtnP();
  virtual ~ZtnP(){}

  int  flavour;//0=mumu,1=ee,2=wrong   float         TnP_pt;   //[nTnP]
  float         TnP_eta;   //[nTnP]
  float         TnP_phi;   //[nTnP]
  float         TnP_mass;   //[nTnP]
  int           TnP_hasFSR;   //[nTnP]
  float         TnP_mll;   //[nTnP]
  int           TnP_l1_pdgId;   //[nTnP]
  float         TnP_l1_pt;   //[nTnP]
  float         TnP_l1_eta;   //[nTnP]
  float         TnP_l1_phi;   //[nTnP]
  float         TnP_l1_mass;   //[nTnP]
  int           TnP_l1_charge;   //[nTnP]
  int           TnP_l1_tightId;   //[nTnP]
  int           TnP_l1_looseId;   //[nTnP]
  float         TnP_l1_dxy;   //[nTnP]
  float         TnP_l1_dz;   //[nTnP]
  float         TnP_l1_edxy;   //[nTnP]
  float         TnP_l1_edz;   //[nTnP]
  float         TnP_l1_ip3d;   //[nTnP]
  float         TnP_l1_sip3d;   //[nTnP]
  float         TnP_l1_ptErr;   //[nTnP]
  int           TnP_l1_lostHits;   //[nTnP]
  int           TnP_l1_trackerLayers;   //[nTnP]
  int           TnP_l1_pixelLayers;   //[nTnP]
  float         TnP_l1_etaSc;   //[nTnP]
  int           TnP_l1_isGap;   //[nTnP]
  float         TnP_l1_r9;   //[nTnP]
  int           TnP_l1_convVeto;   //[nTnP]
  float         TnP_l1_mvaIdSpring15;   //[nTnP]
  float         TnP_l1_relIsoAfterFSR;   //[nTnP]
  float         TnP_l1_chargedHadIso03;   //[nTnP]
  int           TnP_l1_hasOwnFSR;   //[nTnP]
  int           TnP_l1_hlt1L;   //[nTnP]
  float         TnP_l1_p4WithFSR_pt;   //[nTnP]
  float         TnP_l1_p4WithFSR_eta;   //[nTnP]
  float         TnP_l1_p4WithFSR_phi;   //[nTnP]
  float         TnP_l1_p4WithFSR_mass;   //[nTnP]
  int           TnP_l2_pdgId;   //[nTnP]
  float         TnP_l2_pt;   //[nTnP]
  float         TnP_l2_eta;   //[nTnP]
  float         TnP_l2_phi;   //[nTnP]
  float         TnP_l2_mass;   //[nTnP]
  int           TnP_l2_charge;   //[nTnP]
  int           TnP_l2_tightId;   //[nTnP]
  int           TnP_l2_looseId;   //[nTnP]
  float         TnP_l2_dxy;   //[nTnP]
  float         TnP_l2_dz;   //[nTnP]
  float         TnP_l2_edxy;   //[nTnP]
  float         TnP_l2_edz;   //[nTnP]
  float         TnP_l2_ip3d;   //[nTnP]
  float         TnP_l2_sip3d;   //[nTnP]
  float         TnP_l2_ptErr;   //[nTnP]
  int           TnP_l2_lostHits;   //[nTnP]
  int           TnP_l2_trackerLayers;   //[nTnP]
  int           TnP_l2_pixelLayers;   //[nTnP]
  float         TnP_l2_etaSc;   //[nTnP]
  int           TnP_l2_isGap;   //[nTnP]
  float         TnP_l2_r9;   //[nTnP]
  int           TnP_l2_convVeto;   //[nTnP]
  float         TnP_l2_mvaIdSpring15;   //[nTnP]
  float         TnP_l2_relIsoAfterFSR;   //[nTnP]
  float         TnP_l2_chargedHadIso03;   //[nTnP]
  int           TnP_l2_hasOwnFSR;   //[nTnP]
  int           TnP_l2_hlt1L;   //[nTnP]
  float         TnP_l2_p4WithFSR_pt;   //[nTnP]
  float         TnP_l2_p4WithFSR_eta;   //[nTnP]
  float         TnP_l2_p4WithFSR_phi;   //[nTnP]
  float         TnP_l2_p4WithFSR_mass;   //[nTnP]
  ClassDef(ZtnP,1)
};
}
#endif
