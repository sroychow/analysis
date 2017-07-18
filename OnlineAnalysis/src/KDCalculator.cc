#include <iostream>
#include <algorithm>
#include <vector>
#include "AnalysisSpace/OnlineAnalysis/src/KDCalculator.h"
#include "AnalysisSpace/OnlineAnalysis/src/LeptonIsoCalculator.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"

#include <TRandom.h>
KDCalculator::KDCalculator(const bool verbose) {
  verbosity_ = verbose;
}

KDCalculator::~KDCalculator()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

float KDCalculator::computeKDforbestZZ(Mela* mela, vhtm::ZZcandidate& zzcand) {//, std::vector<pat::Jet>& CleanJets) {
  zzcand.setP4Vectors();
  std::vector<TLorentzVector>  lepP4Vec = zzcand.lepP4Vec;
  std::vector<int>  lepids = zzcand.lepcharge;
  int id11 = lepids[0];
  int id12 = lepids[1];
  int id21 = lepids[2];
  int id22 = lepids[3];
  // double ZZMass = zzcand.m4l;

  // Lepton TLorentzVectors, including FSR 
  SimpleParticleCollection_t daughters;
  daughters.push_back(SimpleParticle_t(id11, lepP4Vec[0]));
  daughters.push_back(SimpleParticle_t(id12, lepP4Vec[1]));
  daughters.push_back(SimpleParticle_t(id21, lepP4Vec[2]));
  daughters.push_back(SimpleParticle_t(id22, lepP4Vec[3]));  
  
  SimpleParticleCollection_t associated;
  mela->setInputEvent(&daughters, &associated, 0, 0);
  mela->setCurrentCandidateFromIndex(0);
  float me_0plus_JHU_tmp, me_qqZZ_MCFM_tmp;
  mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  mela->computeP(me_0plus_JHU_tmp, true);            
  mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
  mela->computeP(me_qqZZ_MCFM_tmp, true);
  float D_bkg_kin_tmp = me_0plus_JHU_tmp / (me_0plus_JHU_tmp + me_qqZZ_MCFM_tmp);
  zzcand.Dbkgkin(D_bkg_kin_tmp);
  mela->resetInputEvent(); 
  return D_bkg_kin_tmp;
}
//
// member functions
//
void KDCalculator::computeKDforcategory(Mela* mela, vhtm::ZZcandidate& zzcand,
                                         std::vector<pat::Jet>& cleanJets) {
  std::vector<TLorentzVector>  lepP4VecWithfsr = zzcand.lepP4VecWithfsr;
  std::vector<int>  lepids = zzcand.lepcharge;
  int id11 = lepids[0];
  int id12 = lepids[1];
  int id21 = lepids[2];
  int id22 = lepids[3];
  double mass4l = zzcand.m4l;

  // Lepton TLorentzVectors, including FSR 
  SimpleParticleCollection_t daughters;
  daughters.push_back(SimpleParticle_t(id11, lepP4VecWithfsr[0]));
  daughters.push_back(SimpleParticle_t(id12, lepP4VecWithfsr[1]));
  daughters.push_back(SimpleParticle_t(id21, lepP4VecWithfsr[2]));
  daughters.push_back(SimpleParticle_t(id22, lepP4VecWithfsr[3]));    

  //Add associted particles
  SimpleParticleCollection_t associated; 
  if(cleanJets.size() > 0) associated.push_back(SimpleParticle_t(0, HZZ4lUtil::getP4(cleanJets[0])));
  if(cleanJets.size() > 1) associated.push_back(SimpleParticle_t(0, HZZ4lUtil::getP4(cleanJets[1])));
  mela->setInputEvent(&daughters, &associated, 0, 0);
  mela->setCurrentCandidateFromIndex(0);
  
  // ME
  float D_bkg_kin = 999., D_bkg = 999., D_g4 = 999., D_g1g4 = 999.;
  float me_0plus_JHU=999.0, me_qqZZ_MCFM=999.0, p0plus_m4l=999.0, bkg_m4l=999.0;   
  float p0minus_VAJHU=999.0, pg1g4_VAJHU=999.0, Dgg10_VAMCFM=999.0;
  float phjj_VAJHU=999.0, pvbf_VAJHU=999.0, pAux_vbf_VAJHU=999.0;
  float pwh_hadronic_VAJHU=999.0;
  float pzh_hadronic_VAJHU=999.0;
  float D_HadWH=999.0, D_HadZH=999.0; 
  float D_VBF=999.0, D_VBF1j=999.0;
  float phj_VAJHU;
  //compute all KD
  mela->setInputEvent(&daughters, &associated, 0, 0);
  mela->setCurrentCandidateFromIndex(0);
  
  mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  mela->computeP(me_0plus_JHU, true);
  
  mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
  mela->computeP(p0minus_VAJHU, true);
  
  pg1g4_VAJHU=0.0;
  mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
  (mela->selfDHggcoupl)[0][0][0]=1.;
  (mela->selfDHzzcoupl)[0][0][0]=1.;
  (mela->selfDHzzcoupl)[0][3][0]=1.;
  mela->computeP(pg1g4_VAJHU, true);
  pg1g4_VAJHU -= me_0plus_JHU+p0minus_VAJHU;
  
  mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
  mela->computeP(me_qqZZ_MCFM, true);
  
  mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);
  
  mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);
  
  mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
  mela->computePM4l(TVar::SMSyst_None, bkg_m4l);
  
  D_bkg_kin = me_0plus_JHU/(me_0plus_JHU+me_qqZZ_MCFM*getDbkgkinConstant(id11*id12*id21*id22,mass4l)); 
  D_bkg = me_0plus_JHU*p0plus_m4l/(me_0plus_JHU*p0plus_m4l+me_qqZZ_MCFM*bkg_m4l*getDbkgConstant(id11*id12*id21*id22,mass4l)); // superMELA 
  D_g4 = me_0plus_JHU/(me_0plus_JHU+pow(2.521, 2)*p0minus_VAJHU); 
  D_g1g4 = pg1g4_VAJHU*2.521/(me_0plus_JHU+pow(2.521, 2)*p0minus_VAJHU); 
  zzcand.kdmap["D_bkg_kin"] = D_bkg_kin;
  zzcand.kdmap["D_bkg"] = D_bkg;
  zzcand.kdmap["D_g4"] = D_g4;
  zzcand.kdmap["D_g1g4"] = D_g1g4;

  D_VBF = -1.0; D_HadWH = -1.0; D_HadZH = -1.0;
  if (cleanJets.size() >=2){
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(pvbf_VAJHU, true);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
    mela->computeProdP(phjj_VAJHU, true);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
    mela->computeProdP(pwh_hadronic_VAJHU, true);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
    mela->computeProdP(pzh_hadronic_VAJHU, true);
    D_VBF = pvbf_VAJHU/(pvbf_VAJHU+phjj_VAJHU*getDVBF2jetsConstant(mass4l) ); // VBF(2j) vs. gg->H+2j
    D_HadWH = pwh_hadronic_VAJHU/(pwh_hadronic_VAJHU+1e-3*phjj_VAJHU ); // W(->2j)H vs. gg->H+2j
    D_HadZH = pzh_hadronic_VAJHU/(pzh_hadronic_VAJHU+1e-4*phjj_VAJHU ); // Z(->2j)H vs. gg->H+2j
  } 
  zzcand.kdmap["D_VBF"] = D_VBF;
  zzcand.kdmap["D_HadWH"] = D_HadWH;
  zzcand.kdmap["D_HadZH"] = D_HadZH;

  D_VBF1j = -1;
  if (cleanJets.size() ==1) {
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
    mela->computeProdP(phj_VAJHU, true);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela->computeProdP(pvbf_VAJHU, true); // Un-integrated ME
    mela->getPAux(pAux_vbf_VAJHU); // = Integrated / un-integrated
    D_VBF1j = pvbf_VAJHU*pAux_vbf_VAJHU/(pvbf_VAJHU*pAux_vbf_VAJHU+phj_VAJHU*getDVBF1jetConstant(mass4l)); // VBF(1j) vs. gg->H+1j
  }
  zzcand.kdmap["D_VBF1j"] = D_VBF1j;    

  //KD with Q/J tagger
  float jetPgOverPq[2];
  for (unsigned int j=0; j<cleanJets.size(); j++){
    if(j > 2)  break; //only leading and sub-leading jet
    float QGL = cleanJets[j].userFloat("qgLikelihood");
    if(QGL<0.){ // if the q/g tagger has the error value (-1.), take a random one instead
      TRandom3 rand;
      rand.SetSeed(abs(static_cast<int>(TMath::Sin(cleanJets[j].phi())*100000)));
      QGL = rand.Uniform();
    }
    jetPgOverPq[j] = 1./QGL - 1.;
  }
  float D_VBF1j_QG = -1.; 
  float D_VBF2j_QG = -1.;
  float D_HadWH_QG = -1.;
  float D_HadZH_QG = -1.;
  if(cleanJets.size() == 1) {
    D_VBF1j_QG = 1./(1.+ (1./D_VBF1j - 1.) * pow(jetPgOverPq[0], 1./3.)); 
  } else if(cleanJets.size() >= 2) {
    D_VBF2j_QG = 1./(1.+ (1./D_VBF - 1.) * pow(jetPgOverPq[0]*jetPgOverPq[1], 1./3.));
    D_HadWH_QG = 1./(1.+ (1./D_HadWH - 1.) * pow(jetPgOverPq[0]*jetPgOverPq[1], 1./3.));
    D_HadZH_QG = 1./(1.+ (1./D_HadZH - 1.) * pow(jetPgOverPq[0]*jetPgOverPq[1], 1./3.));
  } 
  zzcand.kdmap["D_VBF1j_QG"] = D_VBF1j_QG;
  zzcand.kdmap["D_VBF2j_QG"] = D_VBF2j_QG;
  zzcand.kdmap["D_HadWH_QG"] = D_HadWH_QG;
  zzcand.kdmap["D_HadZH_QG"] = D_HadZH_QG;
  //dont forget
  mela->resetInputEvent(); 
}
