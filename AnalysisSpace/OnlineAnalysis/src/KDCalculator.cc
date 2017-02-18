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
  /*
  //Dont we need jets in the associated vector?? Difference from Higgs twiki
  //COMMENT:: Confirm with others!!
  SimpleParticleCollection_t associated;
  unsigned int nCandidates = 0;
  for(unsigned int i = 0; i < CleanJets.size(); i ++) {  
    associated.push_back(SimpleParticle_t(0, TLorentzVector(CleanJets[i].pt(), 
                         CleanJets[i].eta(), CleanJets[i].phi(), CleanJets[i].energy())));
    nCandidates++;
  }
    
  mela->setInputEvent(&daughters, &associated, 0, 0);
  mela->setCurrentCandidateFromIndex(0);
    
  float p0plus_VAJHU=0;
  mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
  mela->computeP(p0plus_VAJHU, true);

  float p0minus_VAJHU=0;
  mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
  mela->computeP(p0minus_VAJHU, true);

  float pg1g4_VAJHU=0;
  mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
  (mela->selfDHggcoupl)[0][0][0]=1.;
  (mela->selfDHzzcoupl)[0][0][0]=1.;
  (mela->selfDHzzcoupl)[0][3][0]=1.;
  mela->computeP(pg1g4_VAJHU, true);
  pg1g4_VAJHU -= p0plus_VAJHU+p0minus_VAJHU;

  float bkg_VAMCFM=0;
  mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
  mela->computeP(bkg_VAMCFM, true);
  */

}
//
// member functions
//
