#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TVector2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "AnalysisSpace/OnlineAnalysis/plugins/ZZCandidateProducer.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

//
// constructors and destructor
//
ZZCandidateProducer::ZZCandidateProducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  ZeeTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZeeSrc")),
  ZmumuTag_(iConfig.getUntrackedParameter<edm::InputTag>("ZmumuSrc")),
  ZeeToken_(consumes<std::vector<vhtm::Zee>>(ZeeTag_)),
  ZmumuToken_(consumes<std::vector<vhtm::Zmumu>>(ZmumuTag_)),
  ZZColl_(iConfig.getUntrackedParameter<std::string>("ZZColl", "ZZCandList"))
{
  produces<std::vector<vhtm::ZZcandidate>>(ZZColl_);
}


ZZCandidateProducer::~ZZCandidateProducer()
{
}


// ------------ method called to produce the data  ------------
void
ZZCandidateProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ZZVec_->clear();
  //using namespace edm;
  
  edm::Handle<std::vector<vhtm::Zee>> ZeeVec;
  iEvent.getByToken(ZeeToken_, ZeeVec);

  edm::Handle<std::vector<vhtm::Zmumu>> ZmumuVec;
  iEvent.getByToken(ZmumuToken_, ZmumuVec);
  //4e case  
  if( ZeeVec->size() >= 2 ) {
    for(unsigned int i=0; i<ZeeVec->size();i++) {
      auto& Zi = ZeeVec->at(i);
      for(unsigned int j=i+1; j<ZeeVec->size();j++) {
        auto& Zj = ZeeVec->at(j);
        int whichZ1;
        double mass4l = 99999.;
        if(ZZSelector<vhtm::Zee,vhtm::Zee>(Zi,Zj,whichZ1,mass4l,true) == 0) {
          vhtm::ZZcandidate zztemp;
          zztemp.m4l = mass4l;
          zztemp.flavour = HZZ4lUtil::ZZType::eeee; 
          if(whichZ1 == 1) {
            zztemp.Z1ee = Zi;
            zztemp.Z2ee = Zj;
            zztemp.Z1mdiff = Zi.massdiff;
            zztemp.mZ1 = Zi.mass;
            zztemp.mZ2 = Zj.mass;
            zztemp.Z2lepPtavg = (HZZ4lUtil::getP4(Zj.lep1).Pt() + HZZ4lUtil::getP4(Zj.lep2).Pt())/2.;
          } else {
            zztemp.Z1ee = Zj;
            zztemp.Z2ee = Zi;         
            zztemp.Z1mdiff = Zj.massdiff;
            zztemp.mZ1 = Zj.mass;
            zztemp.mZ2 = Zi.mass;
            zztemp.Z2lepPtavg = (HZZ4lUtil::getP4(Zi.lep1).Pt() + HZZ4lUtil::getP4(Zi.lep2).Pt())/2.;
          } 
          zztemp.setP4Vectors();//very important.used later in KD 
          ZZVec_->push_back(zztemp); 
        }
      }
    }
  }
  
  //4mu case
  if( ZmumuVec->size() >= 2 ) {
    for(unsigned int i=0; i<ZmumuVec->size();i++) {
      auto& Zi = ZmumuVec->at(i);
      for(unsigned int j=i+1; j<ZmumuVec->size();j++) {
        auto& Zj = ZmumuVec->at(j);
        int whichZ1;
        double mass4l = 99999.;
        if(ZZSelector<vhtm::Zmumu,vhtm::Zmumu>(Zi,Zj,whichZ1,mass4l,true) == 0) {
          vhtm::ZZcandidate zztemp;
          zztemp.m4l = mass4l;
          zztemp.flavour = HZZ4lUtil::ZZType::mmmm; 
          if(whichZ1 == 1) {
            zztemp.Z1mm = Zi;
            zztemp.Z2mm = Zj;
            zztemp.Z1mdiff = Zi.massdiff;
            zztemp.mZ1 = Zi.mass;
            zztemp.mZ2 = Zj.mass;
            zztemp.Z2lepPtavg = (HZZ4lUtil::getP4(Zj.lep1).Pt() + HZZ4lUtil::getP4(Zj.lep2).Pt())/2.;
          } else {
            zztemp.Z1mm = Zj;
            zztemp.Z2mm = Zi;         
            zztemp.Z1mdiff = Zj.massdiff;
            zztemp.mZ1 = Zj.mass;
            zztemp.mZ2 = Zi.mass;
            zztemp.Z2lepPtavg = (HZZ4lUtil::getP4(Zi.lep1).Pt() + HZZ4lUtil::getP4(Zi.lep2).Pt())/2.;
          }
          zztemp.setP4Vectors();//very important.used later in KD 
          ZZVec_->push_back(zztemp);  
        }
      }
    }
  }
  
  //2e2mu/2mu2e case
  if( ZmumuVec->size() >= 1 && ZeeVec->size() >=1 ) {
    for(unsigned int i=0; i<ZmumuVec->size();i++) {
      auto& Zi = ZmumuVec->at(i);
      for(unsigned int j=0; j<ZeeVec->size();j++) {
        auto& Zj = ZeeVec->at(j);
        int whichZ1;
        double mass4l = 99999.;
        if(ZZSelector<vhtm::Zmumu,vhtm::Zee>(Zi,Zj,whichZ1,mass4l) == 0) {
          vhtm::ZZcandidate zztemp;
          zztemp.m4l = mass4l;
          if(whichZ1 == 1) {
            zztemp.flavour = HZZ4lUtil::ZZType::mmee; 
            zztemp.Z1mm = Zi;
            zztemp.Z2ee = Zj;
            zztemp.Z1mdiff = Zi.massdiff;
            zztemp.mZ1 = Zi.mass;
            zztemp.mZ2 = Zj.mass;
            zztemp.Z2lepPtavg = (HZZ4lUtil::getP4(Zj.lep1).Pt() + HZZ4lUtil::getP4(Zj.lep2).Pt())/2.;
          } else {
            zztemp.flavour = HZZ4lUtil::ZZType::eemm; 
            zztemp.Z1ee = Zj;
            zztemp.Z2mm = Zi;         
            zztemp.Z1mdiff = Zj.massdiff;
            zztemp.mZ1 = Zj.mass;
            zztemp.mZ2 = Zi.mass;
            zztemp.Z2lepPtavg = (HZZ4lUtil::getP4(Zi.lep1).Pt() + HZZ4lUtil::getP4(Zi.lep2).Pt())/2.;
          } 
          zztemp.setP4Vectors();//very important.used later in KD 
          ZZVec_->push_back(zztemp);  
        }
      }
    }
  }
  // If more than one ZZ candidate survives, choose the one with the Z1 closest in mass to nominal Z. 
  //If two or more candidates include the same Z1 and differ just by the Z2, choose the one with the highest-pT Z2 leptons. 
  //(Other Z pairing choices, e.g. by best Dbkg, could be considered later). 
  if(ZZVec_->size() > 1)
    std::sort(ZZVec_->begin(), ZZVec_->end(), HZZ4lUtil::ZZmasComparator());
  
  std::auto_ptr<std::vector<vhtm::ZZcandidate>> pv1(new std::vector<vhtm::ZZcandidate>(*ZZVec_));
  iEvent.put(pv1,ZZColl_);
}

template <typename T1,typename T2>
int ZZCandidateProducer::ZZSelector(const T1& Zcand1, const T2& Zcand2, int& whichZ1cand,double& m4l, bool sameFlavour) 
{
  //whichZ1cand = 0 for ZZnot selected, =1 if Zcand1=Z1, =2 if Zcand2=Z1
  whichZ1cand = 0;
  // For all possible ZZ pairs, require both Z candidate masses (computed including the FSR photons if present) to be 12 < m(ll(g)) < 120 GeV 
  bool zzm = Zcand1.mass > 12 && Zcand1.mass < 120. && Zcand2.mass > 12 && Zcand2.mass < 120.;
  if(!zzm)  return 1;

  //ΔR(eta,phi)>0.02 between each of the four leptons (to remove ghosts) 
  TLorentzVector z1l1P4 = HZZ4lUtil::getP4(Zcand1.lep1);
  TLorentzVector z1l2P4 = HZZ4lUtil::getP4(Zcand1.lep2);
  TLorentzVector z2l1P4 = HZZ4lUtil::getP4(Zcand2.lep1);
  TLorentzVector z2l2P4 = HZZ4lUtil::getP4(Zcand2.lep2);
  double dra1a2 = z1l1P4.DeltaR(z1l2P4);
  double drb1b2 = z2l1P4.DeltaR(z2l2P4);  
  double dra1b1 = z1l1P4.DeltaR(z2l1P4);
  double dra1b2 = z1l1P4.DeltaR(z2l2P4);
  double dra2b1 = z1l2P4.DeltaR(z2l1P4);
  double dra2b2 = z1l2P4.DeltaR(z2l2P4);
  bool dRlep = dra1a2 > 0.02 && 
               drb1b2 > 0.02 && 
               dra1b1 > 0.02 && 
               dra1b2 > 0.02 && 
               dra2b1 > 0.02 && 
               dra2b2 > 0.02;
  if(!dRlep)    return 2;

  // the two highest-pT leptons of the four pass pT > 20 and 10 GeV 
  std::vector<TLorentzVector> lepP4List;
  lepP4List.push_back(z1l1P4);
  lepP4List.push_back(z1l2P4);
  lepP4List.push_back(z2l1P4);
  lepP4List.push_back(z2l2P4);
  std::sort(lepP4List.begin(), lepP4List.end(), HZZ4lUtil::PtComparatorTL<TLorentzVector>());
  if (lepP4List[0].Pt() <= 20 && lepP4List[1].Pt() <= 10)    return 3;
  
  /*
  //QCD suppression: m(ll)>4 GeV on all four OS pairs that can be made with the four leptons (regardless of flavour). 
  //FSR photons are not used in computing m(ll); since a qcd-induced low mass dilepton (eg jpsi) may have photons nearby 
  //(eg from pi0) and so it's safer to keep the cut on m(ll) alone. 
  */
  TLorentzVector ZaP4, ZbP4,ZafsrP4,ZbfsrP4;
  if (Zcand1.lep1.charge() + Zcand2.lep1.charge() == 0) {
    ZaP4 = z1l1P4 + z2l1P4;
    if(Zcand1.lep1hasfsr) {
      ZafsrP4 += HZZ4lUtil::getP4(Zcand1.fsrl1);
      z1l1P4 +=  HZZ4lUtil::getP4(Zcand1.fsrl1);
    }
    if(Zcand2.lep1hasfsr) {
      ZafsrP4 += HZZ4lUtil::getP4(Zcand2.fsrl1); 
      z2l1P4 +=  HZZ4lUtil::getP4(Zcand2.fsrl1);
    }
    ZbP4 = z1l2P4 + z2l2P4;
    if(Zcand1.lep2hasfsr) {
      ZbfsrP4 += HZZ4lUtil::getP4(Zcand1.fsrl2);
      z1l2P4 +=  HZZ4lUtil::getP4(Zcand1.fsrl2);
    }
    if(Zcand2.lep2hasfsr) {
      ZbfsrP4 += HZZ4lUtil::getP4(Zcand2.fsrl2);
      z2l2P4 +=  HZZ4lUtil::getP4(Zcand2.fsrl2);
    }
  }
  else {
    ZaP4 = z1l1P4 + z2l2P4;
    if(Zcand1.lep1hasfsr) {
      ZafsrP4 += HZZ4lUtil::getP4(Zcand1.fsrl1);
      z1l1P4 +=  HZZ4lUtil::getP4(Zcand1.fsrl1);
    }
    if(Zcand2.lep2hasfsr) {
      ZafsrP4 += HZZ4lUtil::getP4(Zcand2.fsrl2);
      z2l2P4 +=  HZZ4lUtil::getP4(Zcand2.fsrl2);
    }
    ZbP4 = z1l2P4 + z2l1P4;
    if(Zcand1.lep2hasfsr) {
      ZbfsrP4 += HZZ4lUtil::getP4(Zcand1.fsrl2);
      z1l2P4 +=  HZZ4lUtil::getP4(Zcand1.fsrl2);
    }
    if(Zcand2.lep1hasfsr) {
      ZbfsrP4 += HZZ4lUtil::getP4(Zcand2.fsrl1);
      z2l1P4 +=  HZZ4lUtil::getP4(Zcand2.fsrl1);
    }
  }
  if (Zcand1.mass <= 4 || Zcand2.mass <= 4 || ZaP4.M() <= 4. || ZbP4.M() <= 4.)  return 4;
  ZaP4 += ZafsrP4;
  ZbP4 += ZbfsrP4;
  
  //define the Z1 as the one with mass closest to the nominal mZ; require mZ1 > 40 GeV. The other Z is the Z2. 
  bool isZcand1closest = false, isZcand2closest = false;
  double Z1massdiff = 99999.;
  if(Zcand1.massdiff < Zcand2.massdiff) {
    if(Zcand1.mass > 40.) {
      isZcand1closest = true;
      Z1massdiff = Zcand1.massdiff;
    }
  }
  else {
    if(Zcand2.mass > 40.) {
      isZcand2closest = true;
      Z1massdiff = Zcand2.massdiff;
    } 
  }  
  if(!isZcand1closest && !isZcand2closest)     return 5;
  /*
  //additional "smart cut" to reject 4mu/4e pairs where the alternative pairing looks like a on-shell Z+low-mass ll. 
  //(This is required to avoid a background increase when selecting the "best candidate" after all cuts; cf. Simon's slides. 
  //The requirement is: !( |mZa-mZ| < |mZ1-mZ| && mZb<12), where Za and Zb are the mass-sorted alternative pairing Z candidates 
  //(Za being the one closest to the nominal Z mass). FSR photons associated to leptons are included in the mZa/mZb computation. 
  */
  if(sameFlavour) {//Za Zb defined above in QCD supression part
    if (std::fabs(ZbP4.M() - HZZ4lUtil::MZnominal) < std::fabs(ZaP4.M() - HZZ4lUtil::MZnominal)) {
      TLorentzVector lv = ZaP4;
      ZaP4 = ZbP4;
      ZbP4 = lv;
    } 
    if (std::fabs(ZaP4.M() - HZZ4lUtil::MZnominal) < Z1massdiff && ZbP4.M() < 12)     return 6;    
  }
  
  //m(4l) > 70 GeV
  m4l =  (z1l1P4 + z1l2P4 + z2l1P4 + z2l2P4).M();//fsr attached in QCD supression part
  if(m4l <= 70.)  return 7; 
  if(isZcand1closest) whichZ1cand = 1;
  else whichZ1cand = 2;
  return 0;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ZZCandidateProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ZZCandidateProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  ZZCandidateProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  ZZCandidateProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
ZZCandidateProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  ZZVec_= new std::vector<vhtm::ZZcandidate>();
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
ZZCandidateProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZZCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZZCandidateProducer);
