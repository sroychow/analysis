#ifndef __AnalysisSpace_OnlineAnalysis_ZTnpAnalyzer_h
#define __AnalysisSpace_OnlineAnalysis_ZTnpAnalyzer_h

/*
\class ZTnpAnalyzer ZTnpAnalyzer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Suvankar Roy Chowdhury
//         Created:  Mon, 06 Feb 2017 11:56:27 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ZTnpAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZTnpAnalyzer(const edm::ParameterSet&);
      ~ZTnpAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      bool verbosity_;
      bool isMC_;
      const edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
      const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      //Selected Electrons token
      const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      //variables for ZtoMuMu TnP tree
      std::string outFileName_;
      TTree* outTree_;
      TFile* outTreeFile_;
      const static Int_t kMaxTnP = 10;
  //////////////////////////////////////////////////////////////
   int nTnP;
   float         TnP_pt[kMaxTnP];   
   float         TnP_eta[kMaxTnP];   
   float         TnP_phi[kMaxTnP];   
   float         TnP_mass[kMaxTnP];   
   int           TnP_hasFSR[kMaxTnP];   
   float         TnP_mll[kMaxTnP];   
   int           TnP_l1_pdgId[kMaxTnP];   
   float         TnP_l1_pt[kMaxTnP];   
   float         TnP_l1_eta[kMaxTnP];   
   float         TnP_l1_phi[kMaxTnP];   
   float         TnP_l1_mass[kMaxTnP];   
   int           TnP_l1_charge[kMaxTnP];   
   int           TnP_l1_tightId[kMaxTnP];   
   int           TnP_l1_looseId[kMaxTnP];   
   float         TnP_l1_dxy[kMaxTnP];   
   float         TnP_l1_dz[kMaxTnP];   
   float         TnP_l1_edxy[kMaxTnP];   
   float         TnP_l1_edz[kMaxTnP];   
   float         TnP_l1_ip3d[kMaxTnP];   
   float         TnP_l1_sip3d[kMaxTnP];   
   float         TnP_l1_ptErr[kMaxTnP];   
   int           TnP_l1_lostHits[kMaxTnP];   
   int           TnP_l1_trackerLayers[kMaxTnP];   
   int           TnP_l1_pixelLayers[kMaxTnP];   
   float         TnP_l1_etaSc[kMaxTnP];   
   int           TnP_l1_isGap[kMaxTnP];   
   float         TnP_l1_r9[kMaxTnP];   
   int           TnP_l1_convVeto[kMaxTnP];   
   float         TnP_l1_mvaIdSpring15[kMaxTnP];   
   float         TnP_l1_relIsoAfterFSR[kMaxTnP];   
   float         TnP_l1_chargedHadIso03[kMaxTnP];   
   int           TnP_l1_hasOwnFSR[kMaxTnP]; 

   int           TnP_l1_mcMatchId[kMaxTnP];   //[nTnP]
   int           TnP_l1_mcMatchAny[kMaxTnP];   //[nTnP]
   float         TnP_l1_mcPt[kMaxTnP];   //[nTnP]
   float         TnP_l1_mcPt1[kMaxTnP]; 
 
   int           TnP_l1_hlt1L[kMaxTnP];   
   float         TnP_l1_p4WithFSR_pt[kMaxTnP];   
   float         TnP_l1_p4WithFSR_eta[kMaxTnP];   
   float         TnP_l1_p4WithFSR_phi[kMaxTnP];   
   float         TnP_l1_p4WithFSR_mass[kMaxTnP];   
   int           TnP_l2_pdgId[kMaxTnP];   
   float         TnP_l2_pt[kMaxTnP];   
   float         TnP_l2_eta[kMaxTnP];   
   float         TnP_l2_phi[kMaxTnP];   
   float         TnP_l2_mass[kMaxTnP];   
   int           TnP_l2_charge[kMaxTnP];   
   int           TnP_l2_tightId[kMaxTnP];   
   int           TnP_l2_looseId[kMaxTnP];   
   float         TnP_l2_dxy[kMaxTnP];   
   float         TnP_l2_dz[kMaxTnP];   
   float         TnP_l2_edxy[kMaxTnP];   
   float         TnP_l2_edz[kMaxTnP];   
   float         TnP_l2_ip3d[kMaxTnP];   
   float         TnP_l2_sip3d[kMaxTnP];   
   float         TnP_l2_ptErr[kMaxTnP];   
   int           TnP_l2_lostHits[kMaxTnP];   
   int           TnP_l2_trackerLayers[kMaxTnP];   
   int           TnP_l2_pixelLayers[kMaxTnP];   
   float         TnP_l2_etaSc[kMaxTnP];   
   int           TnP_l2_isGap[kMaxTnP];   
   float         TnP_l2_r9[kMaxTnP];   
   int           TnP_l2_convVeto[kMaxTnP];   
   float         TnP_l2_mvaIdSpring15[kMaxTnP];   
   float         TnP_l2_relIsoAfterFSR[kMaxTnP];   
   float         TnP_l2_chargedHadIso03[kMaxTnP];   
   int           TnP_l2_hasOwnFSR[kMaxTnP];

   int           TnP_l2_mcMatchId[kMaxTnP];   //[nTnP]
   int           TnP_l2_mcMatchAny[kMaxTnP];   //[nTnP]
   float         TnP_l2_mcPt[kMaxTnP];   //[nTnP]
   float         TnP_l2_mcPt1[kMaxTnP];   

   int           TnP_l2_hlt1L[kMaxTnP];   
   float         TnP_l2_p4WithFSR_pt[kMaxTnP];   
   float         TnP_l2_p4WithFSR_eta[kMaxTnP];   
   float         TnP_l2_p4WithFSR_phi[kMaxTnP];   
   float         TnP_l2_p4WithFSR_mass[kMaxTnP];
  //////////////////////////////////////////////////////////////
      

};
#endif
