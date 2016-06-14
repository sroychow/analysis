#include <memory>
#include<iomanip>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "AnalysisSpace/OnlineAnalysis/plugins/ZTnpAnalyzer.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"


//header file for tree
#include "AnalysisSpace/TreeMaker/interface/Utility.h"



//
// constructors and destructor
//
ZTnpAnalyzer::ZTnpAnalyzer(const edm::ParameterSet& iConfig) :
  fsrTag_(iConfig.getUntrackedParameter<edm::InputTag>("fsrCol",edm::InputTag("packedPFCandidates"))), 
  looseSIPEleFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("looseSIPElectronFSRSrc")),
  looseSIPMuFSRpairTag_(iConfig.getUntrackedParameter<edm::InputTag>("looseSIPMuonFSRSrc")),
  fixedGridRhoFastjetAllTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllTag",edm::InputTag("fixedGridRhoFastjetAll"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  fsrToken_(consumes<pat::PackedCandidateCollection>(fsrTag_)),
  looseElectronFSRToken_(consumes<EleFSRCollection>(looseSIPEleFSRpairTag_)),
  looseMuonFSRToken_(consumes<MuFSRCollection>(looseSIPMuFSRpairTag_)),
  fixedGridRhoFastjetAllToken_(consumes<double>(fixedGridRhoFastjetAllTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  isData_(iConfig.getUntrackedParameter<bool>("isData",false))
{
 
}


ZTnpAnalyzer::~ZTnpAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZTnpAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  zlltnp_->clear();
  //get loose + SIP Leptons + FSR collection 
  edm::Handle<pat::PackedCandidateCollection> fsr;
  iEvent.getByToken(fsrToken_,fsr);
  edm::Handle<EleFSRCollection> looseSIPEleFSRpair;
  iEvent.getByToken(looseElectronFSRToken_,looseSIPEleFSRpair);
  edm::Handle<MuFSRCollection> looseSIPMuFSRpair;
  iEvent.getByToken(looseMuonFSRToken_,looseSIPMuFSRpair);
  
  edm::Handle<double> fixedGridRhoFastjetAll;
  iEvent.getByToken(fixedGridRhoFastjetAllToken_,fixedGridRhoFastjetAll);
  const double evrho = *fixedGridRhoFastjetAll;

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vertexToken_, primaryVertices);
  
  //std::cout << "ZTnP::"
  //         << "PVvalid=" << primaryVertices.isValid()
  //         << "\t\n#Mu=" << looseSIPMuFSRpair->size()
  //         << "\t\n#Ele=" << looseSIPEleFSRpair->size()
  //         << std::endl;
  if (primaryVertices.isValid()) {
    const reco::Vertex& vit = primaryVertices->front(); 
    fillZeeTnp(looseSIPEleFSRpair,fsr,vit,evrho);
    fillZmumuTnp(looseSIPMuFSRpair,fsr,vit);
  }
  std::cout << "ZTnP::ZTnP Candidates=" << zlltnp_->size() << std::endl;
  

}

//select Zmuon TnP candidate
//tag = muon with tightID + ISO + SingleMu ISO 20Trigger
void ZTnpAnalyzer::fillZmumuTnp(const edm::Handle<MuFSRCollection>& looseSIPMuFSRpair, const edm::Handle<pat::PackedCandidateCollection>& fsr,const reco::Vertex& vit) {
  for(unsigned int ti = 0; ti < looseSIPMuFSRpair->size(); ti++) {
    const auto& tag = looseSIPMuFSRpair->at(ti);
    //tag tight mu cut
    if(!tag.first.isPFMuon())    continue; 
    double tagiso = HZZ4lUtil::computeMuonReliso(tag.first, fsr);
    //tag iso cut
    if( tagiso >= 0.35)      continue;
    TLorentzVector tagP4 = HZZ4lUtil::getP4(tag.first);
    TLorentzVector tagP4wfsr = tagP4;
    bool taghasfsr = false;
    if(!tag.second.empty()) {
      tagP4wfsr += HZZ4lUtil::getP4(tag.second.at(0));
      taghasfsr = true;
    }
    for(unsigned int pi = ti+1; pi < looseSIPMuFSRpair->size(); pi++) {
      const auto& probe = looseSIPMuFSRpair->at(pi);
      //opposite charge cut
      if(tag.first.charge() + probe.first.charge() != 0)    continue;
      //probe tight cut
      if(!probe.first.isPFMuon())    continue; 
      TLorentzVector probeP4 = HZZ4lUtil::getP4(probe.first);
      TLorentzVector probeP4wfsr = probeP4;
      bool probehasfsr = false;
      if(!probe.second.empty()) {
        probeP4wfsr += HZZ4lUtil::getP4(probe.second.at(0));
        probehasfsr = true;
      }
      double probeiso = HZZ4lUtil::computeMuonReliso(probe.first, fsr);

      TLorentzVector tnpP4 = tagP4 + probeP4;
      TLorentzVector tnpP4wfsr = tagP4wfsr + probeP4wfsr;

      vhtm::ZtnP  ztemp;
      ztemp.flavour = HZZ4lUtil::ZType::mumu;
      ztemp.TnP_eta = tnpP4.Eta();   
      ztemp.TnP_phi = tnpP4.Phi();   
      ztemp.TnP_mass = tnpP4.M();   //without fsr?
      ztemp.TnP_hasFSR = taghasfsr || probehasfsr;   
      ztemp.TnP_mll = tnpP4wfsr.M();   //with fsr?
      ztemp.TnP_l1_pdgId = (tag.first.charge() > 0) ? 13 : -13;   
      ztemp.TnP_l1_pt = tagP4.Pt();   
      ztemp.TnP_l1_eta = tagP4.Eta();   
      ztemp.TnP_l1_phi = tagP4.Phi();   
      ztemp.TnP_l1_mass = tagP4.M();   
      ztemp.TnP_l1_charge = tag.first.charge();   
      ztemp.TnP_l1_tightId = 1;   
      ztemp.TnP_l1_looseId = 1;   
      ztemp.TnP_l1_dxy = tag.first.muonBestTrack()->dxy(vit.position());   
      ztemp.TnP_l1_dz  = tag.first.muonBestTrack()->dxy(vit.position());;   
      ztemp.TnP_l1_edxy = tag.first.muonBestTrack()->dxyError();   //correct? 
      ztemp.TnP_l1_edz = tag.first.muonBestTrack()->dzError();
      ztemp.TnP_l1_ip3d = tag.first.dB(pat::Muon::PV3D);   //check
      ztemp. TnP_l1_sip3d = tag.first.edB(pat::Muon::PV3D);   //check
      ztemp.TnP_l1_ptErr = tag.first.muonBestTrack()->ptError();   
      ztemp.TnP_l1_lostHits = tag.first.muonBestTrack()->numberOfLostHits();   
      ztemp.TnP_l1_trackerLayers = tag.first.muonBestTrack()->hitPattern().numberOfValidTrackerHits();   
      ztemp.TnP_l1_pixelLayers = tag.first.muonBestTrack()->hitPattern().numberOfValidPixelHits();;   
      //ztemp.TnP_l1_etaSc = -1;   
      //ztemp.TnP_l1_isGap = -1;   
      //ztemp.TnP_l1_r9 = -1;   
      //ztemp.TnP_l1_convVeto = -1;   
      //ztemp.TnP_l1_mvaIdSpring15 = -1;;   
      ztemp.TnP_l1_relIsoAfterFSR = tagiso;   
      ztemp.TnP_l1_chargedHadIso03 = tag.first.pfIsolationR03().sumChargedHadronPt;   
      ztemp.TnP_l1_hasOwnFSR  = 1; //doubt
      ztemp.TnP_l1_hlt1L = 1; //doubt  
      ztemp.TnP_l1_p4WithFSR_pt = tagP4wfsr.Pt();   
      ztemp.TnP_l1_p4WithFSR_eta = tagP4wfsr.Eta();   
      ztemp.TnP_l1_p4WithFSR_phi = tagP4wfsr.Phi();   
      ztemp.TnP_l1_p4WithFSR_mass = tagP4wfsr.M();   

      ztemp.TnP_l2_pdgId = (probe.first.charge() > 0) ? 13 : -13;   
      ztemp.TnP_l2_pt = probeP4.Pt();   
      ztemp.TnP_l2_eta = probeP4.Eta();   
      ztemp.TnP_l2_phi = probeP4.Phi();   
      ztemp.TnP_l2_mass = probeP4.M();   
      ztemp.TnP_l2_charge = probe.first.charge();   
      ztemp.TnP_l2_tightId = 1;   
      ztemp.TnP_l2_looseId = 1;   
      ztemp.TnP_l2_dxy = probe.first.muonBestTrack()->dxy(vit.position());   
      ztemp.TnP_l2_dz  = probe.first.muonBestTrack()->dxy(vit.position());;   
      ztemp.TnP_l2_edxy = probe.first.muonBestTrack()->dxyError();   //correct? 
      ztemp.TnP_l2_edz = probe.first.muonBestTrack()->dzError();
      ztemp.TnP_l2_ip3d = probe.first.dB(pat::Muon::PV3D);   //check
      ztemp. TnP_l2_sip3d = probe.first.edB(pat::Muon::PV3D);   //check
      ztemp.TnP_l2_ptErr = probe.first.muonBestTrack()->ptError();   
      ztemp.TnP_l2_lostHits = probe.first.muonBestTrack()->numberOfLostHits();   
      ztemp.TnP_l2_trackerLayers = probe.first.muonBestTrack()->hitPattern().numberOfValidTrackerHits();   
      ztemp.TnP_l2_pixelLayers = probe.first.muonBestTrack()->hitPattern().numberOfValidPixelHits();;   
      //ztemp.TnP_l2_etaSc = -1;   
      //ztemp.TnP_l2_isGap = -1;   
      //ztemp.TnP_l2_r9 = -1;   
      //ztemp.TnP_l2_convVeto = -1;   
      //ztemp.TnP_l2_mvaIdSpring15 = -1;;   
      ztemp.TnP_l2_relIsoAfterFSR = probeiso;   
      ztemp.TnP_l2_chargedHadIso03 = probe.first.pfIsolationR03().sumChargedHadronPt;   
      ztemp.TnP_l2_hasOwnFSR  = 1; //doubt
      ztemp.TnP_l2_hlt1L = 1; //doubt  
      ztemp.TnP_l2_p4WithFSR_pt = probeP4wfsr.Pt();   
      ztemp.TnP_l2_p4WithFSR_eta = probeP4wfsr.Eta();   
      ztemp.TnP_l2_p4WithFSR_phi = probeP4wfsr.Phi();   
      ztemp.TnP_l2_p4WithFSR_mass = probeP4wfsr.M();  
      zlltnp_->push_back(ztemp); 
    }
  }
}

void 
ZTnpAnalyzer::fillZeeTnp(const edm::Handle<EleFSRCollection>& looseSIPEleFSRpair, const edm::Handle<pat::PackedCandidateCollection>& fsr,
                         const reco::Vertex& vit, const double& evrho) 
{
  for(unsigned int ti = 0; ti < looseSIPEleFSRpair->size(); ti++) {
    const auto& tag = looseSIPEleFSRpair->at(ti);
    //tag tight ele cut
    if(!HZZ4lUtil::passBDT(std::fabs(tag.first.superCluster()->eta()), tag.first.pt(), 
                  tag.first.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")))    continue; 
    double tagiso = HZZ4lUtil::computeElectronReliso(tag.first, fsr,evrho);
    //tag iso cut
    if( tagiso >= 0.35)      continue;
    TLorentzVector tagP4 = HZZ4lUtil::getP4(tag.first);
    TLorentzVector tagP4wfsr = tagP4;
    bool taghasfsr = false;
    if(!tag.second.empty()) {
      tagP4wfsr += HZZ4lUtil::getP4(tag.second.at(0));
      taghasfsr = true;
    }
    for(unsigned int pi = ti+1; pi < looseSIPEleFSRpair->size(); pi++) {
      const auto& probe = looseSIPEleFSRpair->at(pi);
      //opposite charge cut
      if(tag.first.charge() + probe.first.charge() != 0)    continue;
      //probe tight cut
      if(!!HZZ4lUtil::passBDT(std::fabs(probe.first.superCluster()->eta()), probe.first.pt(), 
                  probe.first.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")))    continue; 
      TLorentzVector probeP4 = HZZ4lUtil::getP4(probe.first);
      TLorentzVector probeP4wfsr = probeP4;
      bool probehasfsr = false;
      if(!probe.second.empty()) {
        probeP4wfsr += HZZ4lUtil::getP4(probe.second.at(0));
        probehasfsr = true;
      }
      double probeiso = HZZ4lUtil::computeElectronReliso(probe.first, fsr,evrho);

      TLorentzVector tnpP4 = tagP4 + probeP4;
      TLorentzVector tnpP4wfsr = tagP4wfsr + probeP4wfsr;

      vhtm::ZtnP  ztemp;
      ztemp.flavour = HZZ4lUtil::ZType::ee;
      ztemp.TnP_eta = tnpP4.Eta();   
      ztemp.TnP_phi = tnpP4.Phi();   
      ztemp.TnP_mass = tnpP4.M();   //without fsr?
      ztemp.TnP_hasFSR = taghasfsr || probehasfsr;   
      ztemp.TnP_mll = tnpP4wfsr.M();   //with fsr?
      ztemp.TnP_l1_pdgId = (tag.first.charge() > 0) ? 11 : -11;   
      ztemp.TnP_l1_pt = tagP4.Pt();   
      ztemp.TnP_l1_eta = tagP4.Eta();   
      ztemp.TnP_l1_phi = tagP4.Phi();   
      ztemp.TnP_l1_mass = tagP4.M();   
      ztemp.TnP_l1_charge = tag.first.charge();   
      ztemp.TnP_l1_tightId = 1;   
      ztemp.TnP_l1_looseId = 1;   
      ztemp.TnP_l1_dxy = tag.first.gsfTrack()->dxy(vit.position());   
      ztemp.TnP_l1_dz  = tag.first.gsfTrack()->dxy(vit.position());;   
      ztemp.TnP_l1_edxy = tag.first.gsfTrack()->dxyError();   //correct? 
      ztemp.TnP_l1_edz = tag.first.gsfTrack()->dzError();
      ztemp.TnP_l1_ip3d = tag.first.dB(pat::Electron::PV3D);   //check
      ztemp. TnP_l1_sip3d = tag.first.edB(pat::Electron::PV3D);   //check
      ztemp.TnP_l1_ptErr = tag.first.gsfTrack()->ptError();   
      ztemp.TnP_l1_lostHits = tag.first.gsfTrack()->numberOfLostHits();   
      ztemp.TnP_l1_trackerLayers = tag.first.gsfTrack()->hitPattern().numberOfValidTrackerHits();   
      ztemp.TnP_l1_pixelLayers = tag.first.gsfTrack()->hitPattern().numberOfValidPixelHits();;   
      ztemp.TnP_l1_etaSc = tag.first.superCluster()->eta();   
      ztemp.TnP_l1_isGap = tag.first.isGap();   
      ztemp.TnP_l1_r9 = -1;   
      ztemp.TnP_l1_convVeto = tag.first.passConversionVeto();   
      ztemp.TnP_l1_mvaIdSpring15 = tag.first.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");   
      ztemp.TnP_l1_relIsoAfterFSR = tagiso;   
      ztemp.TnP_l1_chargedHadIso03 = tag.first.pfIsolationVariables().sumChargedHadronPt;   
      ztemp.TnP_l1_hasOwnFSR  = 1; //doubt
      ztemp.TnP_l1_hlt1L = 1; //doubt  
      ztemp.TnP_l1_p4WithFSR_pt = tagP4wfsr.Pt();   
      ztemp.TnP_l1_p4WithFSR_eta = tagP4wfsr.Eta();   
      ztemp.TnP_l1_p4WithFSR_phi = tagP4wfsr.Phi();   
      ztemp.TnP_l1_p4WithFSR_mass = tagP4wfsr.M();   

      ztemp.TnP_l2_pdgId = (probe.first.charge() > 0) ? 11 : -11;   
      ztemp.TnP_l2_pt = probeP4.Pt();   
      ztemp.TnP_l2_eta = probeP4.Eta();   
      ztemp.TnP_l2_phi = probeP4.Phi();   
      ztemp.TnP_l2_mass = probeP4.M();   
      ztemp.TnP_l2_charge = probe.first.charge();   
      ztemp.TnP_l2_tightId = 1;   
      ztemp.TnP_l2_looseId = 1;   
      ztemp.TnP_l2_dxy = probe.first.gsfTrack()->dxy(vit.position());   
      ztemp.TnP_l2_dz  = probe.first.gsfTrack()->dxy(vit.position());;   
      ztemp.TnP_l2_edxy = probe.first.gsfTrack()->dxyError();   //correct? 
      ztemp.TnP_l2_edz = probe.first.gsfTrack()->dzError();
      ztemp.TnP_l2_ip3d = probe.first.dB(pat::Electron::PV3D);   //check
      ztemp. TnP_l2_sip3d = probe.first.edB(pat::Electron::PV3D);   //check sip or error on sip?
      ztemp.TnP_l2_ptErr = probe.first.gsfTrack()->ptError();   
      ztemp.TnP_l2_lostHits = probe.first.gsfTrack()->numberOfLostHits();   
      ztemp.TnP_l2_trackerLayers = probe.first.gsfTrack()->hitPattern().numberOfValidTrackerHits();   
      ztemp.TnP_l2_pixelLayers = probe.first.gsfTrack()->hitPattern().numberOfValidPixelHits();;   
      ztemp.TnP_l2_etaSc = probe.first.superCluster()->eta();   
      ztemp.TnP_l2_isGap = probe.first.isGap();   
      ztemp.TnP_l2_r9 = -1;   
      ztemp.TnP_l2_convVeto = probe.first.passConversionVeto();   
      ztemp.TnP_l2_mvaIdSpring15 = probe.first.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");   
      ztemp.TnP_l2_relIsoAfterFSR = probeiso;   
      ztemp.TnP_l2_chargedHadIso03 = probe.first.pfIsolationVariables().sumChargedHadronPt;   
      ztemp.TnP_l2_hasOwnFSR  = 1; //doubt
      ztemp.TnP_l2_hlt1L = 1; //doubt  
      ztemp.TnP_l2_p4WithFSR_pt = probeP4wfsr.Pt();   
      ztemp.TnP_l2_p4WithFSR_eta = probeP4wfsr.Eta();   
      ztemp.TnP_l2_p4WithFSR_phi = probeP4wfsr.Phi();   
      ztemp.TnP_l2_p4WithFSR_mass = probeP4wfsr.M();  
      zlltnp_->push_back(ztemp); 
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZTnpAnalyzer::beginJob()
{
  zlltnp_ = new std::vector<vhtm::ZtnP>();
  TTree* tree = vhtm::Utility::getTree("vhtree");
  tree->Branch("ZTnPCand", "std::vector<vhtm::ZtnP>", &zlltnp_, 32000, -1);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZTnpAnalyzer::endJob() 
{

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZTnpAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZTnpAnalyzer);
