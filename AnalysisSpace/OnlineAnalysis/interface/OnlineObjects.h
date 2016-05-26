#ifndef __AnalysisSpace_OnlineAnalysis_OnlineObjects_h
#define __AnalysisSpace_OnlineAnalysis_OnlineObjects_h

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "TLorentzVector.h"

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
class SelectedEvent {
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
    void reset();
    ClassDef(SelectedEvent,1)
};
}
#endif
