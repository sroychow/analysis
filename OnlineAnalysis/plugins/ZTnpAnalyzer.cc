#include <iostream>
#include <algorithm>
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "AnalysisSpace/OnlineAnalysis/plugins/ZTnpAnalyzer.h"
#include "AnalysisSpace/OnlineAnalysis/src/PhysicsObjectSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"
#include "AnalysisSpace/OnlineAnalysis/src/LeptonIsoCalculator.h"
// constructors and destructor
//
ZTnpAnalyzer::ZTnpAnalyzer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  syncFilename_(iConfig.getUntrackedParameter<std::string>("syncFilename")),
  fixedGridRhoFastjetAllToken_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllTag",edm::InputTag("fixedGridRhoFastjetAll")))),  
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcandSrc"))),
  outFileName_((iConfig.getParameter<std::string>("tnpFile"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
  edm::Service<TFileService> outFile;
  outTree_ = outFile->make<TTree>("zmmtree","muTnP tree for iso");
   //Set TnP tree 
  outTreeFile_->cd();
  outTree_ = new TTree("tree","muTnP tree for iso");
  outTree_->Branch("run",&run,"run/I");
  outTree_->Branch("lumi",&lumi,"lumi/I");
  outTree_->Branch("event",&event,"event/I");
  outTree_->Branch("isData",&isData,"isData/I");
  outTree_->Branch("nTnP",&nTnP,"nTnP/I");
  
  outTree_->Branch("TnP_pt",TnP_pt,"TnP_pt[nTnP]/F");
  outTree_->Branch("TnP_eta",TnP_eta,"TnP_eta[nTnP]/F");
  outTree_->Branch("TnP_phi",TnP_phi,"TnP_phi[nTnP]/F");
  outTree_->Branch("TnP_mass",TnP_mass,"TnP_mass[nTnP]/F");
  outTree_->Branch("TnP_hasFSR",TnP_hasFSR,"TnP_hasFSR[nTnP]/I");
  outTree_->Branch("TnP_mll",TnP_mll,"TnP_mll[nTnP]/F");
  outTree_->Branch("TnP_l1_pdgId",TnP_l1_pdgId,"TnP_l1_pdgId[nTnP]/I");
  outTree_->Branch("TnP_l1_pt",TnP_l1_pt,"TnP_l1_pt[nTnP]/F");
  outTree_->Branch("TnP_l1_eta",TnP_l1_eta,"TnP_l1_eta[nTnP]/F");
  outTree_->Branch("TnP_l1_phi",TnP_l1_phi,"TnP_l1_phi[nTnP]/F");
  outTree_->Branch("TnP_l1_mass",TnP_l1_mass,"TnP_l1_mass[nTnP]/F");
  outTree_->Branch("TnP_l1_charge",TnP_l1_charge,"TnP_l1_charge[nTnP]/I");
  outTree_->Branch("TnP_l1_tightId",TnP_l1_tightId,"TnP_l1_tightId[nTnP]/I");
  outTree_->Branch("TnP_l1_looseId",TnP_l1_looseId,"TnP_l1_looseId[nTnP]/I");
  outTree_->Branch("TnP_l1_dxy",TnP_l1_dxy,"TnP_l1_dxy[nTnP]/F");
  outTree_->Branch("TnP_l1_dz",TnP_l1_dz,"TnP_l1_dz[nTnP]/F");
  outTree_->Branch("TnP_l1_edxy",TnP_l1_edxy,"TnP_l1_edxy[nTnP]/F");
  outTree_->Branch("TnP_l1_edz",TnP_l1_edz,"TnP_l1_edz[nTnP]/F");
  outTree_->Branch("TnP_l1_ip3d",TnP_l1_ip3d,"TnP_l1_ip3d[nTnP]/F");
  outTree_->Branch("TnP_l1_sip3d",TnP_l1_sip3d,"TnP_l1_sip3d[nTnP]/F");
  outTree_->Branch("TnP_l1_ptErr",TnP_l1_ptErr,"TnP_l1_ptErr[nTnP]/F");
  outTree_->Branch("TnP_l1_lostHits",TnP_l1_lostHits,"TnP_l1_lostHits[nTnP]/I");
  outTree_->Branch("TnP_l1_trackerLayers",TnP_l1_trackerLayers,"TnP_l1_trackerLayer[nTnP]/I");
  outTree_->Branch("TnP_l1_pixelLayers",TnP_l1_pixelLayers,"TnP_l1_pixelLayers[nTnP]/I");
  outTree_->Branch("TnP_l1_etaSc",TnP_l1_etaSc,"TnP_l1_etaSc[nTnP]/F");
  outTree_->Branch("TnP_l1_isGap",TnP_l1_isGap,"TnP_l1_isGap[nTnP]/F");
  outTree_->Branch("TnP_l1_r9",TnP_l1_r9,"TnP_l1_r9[nTnP]/F");
  outTree_->Branch("TnP_l1_convVeto",TnP_l1_convVeto,"TnP_l1_convVeto[nTnP]/F");
  outTree_->Branch("TnP_l1_mvaIdSpring15",TnP_l1_mvaIdSpring15,"TnP_l1_mvaIdSpring15[nTnP]/F");
  outTree_->Branch("TnP_l1_relIsoAfterFSR",TnP_l1_relIsoAfterFSR,"TnP_l1_relIsoAfterFSR[nTnP]/F");
  outTree_->Branch("TnP_l1_chargedHadIso03",TnP_l1_chargedHadIso03,"TnP_l1_chargedHadIso03[nTnP]/F");
  outTree_->Branch("TnP_l1_hasOwnFSR",TnP_l1_hasOwnFSR,"TnP_l1_hasOwnFSR[nTnP]/I");
  if(isMC()) {
    outTree_->Branch("TnP_l1_mcMatchId",TnP_l1_mcMatchId,"TnP_l1_mcMatchId[nTnP]/I");
    outTree_->Branch("TnP_l1_mcMatchAny",TnP_l1_mcMatchAny,"TnP_l1_mcMatchAny[nTnP]/I");
    outTree_->Branch("TnP_l1_mcPt",TnP_l1_mcPt,"TnP_l1_mcPt[nTnP]/F");
    outTree_->Branch("TnP_l1_mcPt1",TnP_l1_mcPt1,"TnP_l1_mcPt1[nTnP]/F");
  }

  outTree_->Branch("TnP_l1_hlt1L",TnP_l1_hlt1L,"TnP_l1_hlt1L[nTnP]/I");
  outTree_->Branch("TnP_l1_p4WithFSR_pt",TnP_l1_p4WithFSR_pt,"TnP_l1_p4WithFSR_pt[nTnP]/F");
  outTree_->Branch("TnP_l1_p4WithFSR_eta",TnP_l1_p4WithFSR_eta,"TnP_l1_p4WithFSR_eta[nTnP]/F");
  outTree_->Branch("TnP_l1_p4WithFSR_phi",TnP_l1_p4WithFSR_phi,"TnP_l1_p4WithFSR_phi[nTnP]/F");
  outTree_->Branch("TnP_l1_p4WithFSR_mass",TnP_l1_p4WithFSR_mass,"TnP_l1_p4WithFSR_mass[nTnP]/F");

  outTree_->Branch("TnP_l2_pdgId",TnP_l2_pdgId,"TnP_l2_pdgId[nTnP]/I");
  outTree_->Branch("TnP_l2_pt",TnP_l2_pt,"TnP_l2_pt[nTnP]/F");
  outTree_->Branch("TnP_l2_eta",TnP_l2_eta,"TnP_l2_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_phi",TnP_l2_phi,"TnP_l2_phi[nTnP]/F");
  outTree_->Branch("TnP_l2_mass",TnP_l2_mass,"TnP_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_charge",TnP_l2_charge,"TnP_l2_charge[nTnP]/I");
  outTree_->Branch("TnP_l2_tightId",TnP_l2_tightId,"TnP_l2_tightId[nTnP]/I");
  outTree_->Branch("TnP_l2_looseId",TnP_l2_looseId,"TnP_l2_looseId[nTnP]/I");
  outTree_->Branch("TnP_l2_dxy",TnP_l2_dxy,"TnP_l2_dxy[nTnP]/F");
  outTree_->Branch("TnP_l2_dz",TnP_l2_dz,"TnP_l2_dz[nTnP]/F");
  outTree_->Branch("TnP_l2_edxy",TnP_l2_edxy,"TnP_l2_edxy[nTnP]/F");
  outTree_->Branch("TnP_l2_edz",TnP_l2_edz,"TnP_l2_edz[nTnP]/F");
  outTree_->Branch("TnP_l2_ip3d",TnP_l2_ip3d,"TnP_l2_ip3d[nTnP]/F");
  outTree_->Branch("TnP_l2_sip3d",TnP_l2_sip3d,"TnP_l2_sip3d[nTnP]/F");
  outTree_->Branch("TnP_l2_ptErr",TnP_l2_ptErr,"TnP_l2_ptErr[nTnP]/F");
  outTree_->Branch("TnP_l2_lostHits",TnP_l2_lostHits,"TnP_l2_lostHits[nTnP]/I");
  outTree_->Branch("TnP_l2_trackerLayers",TnP_l2_trackerLayers,"TnP_l2_trackerLayers[nTnP]/I");
  outTree_->Branch("TnP_l2_pixelLayers",TnP_l2_pixelLayers,"TnP_l2_pixelLayers[nTnP]/I");
  outTree_->Branch("TnP_l2_etaSc",TnP_l2_etaSc,"TnP_l2_etaSc[nTnP]/F");
  outTree_->Branch("TnP_l2_isGap",TnP_l2_isGap,"TnP_l2_isGap[nTnP]/I");
  outTree_->Branch("TnP_l2_r9",TnP_l2_r9,"TnP_l2_r9[nTnP]/F");
  outTree_->Branch("TnP_l2_convVeto",TnP_l2_convVeto,"TnP_l2_convVeto[nTnP]/I");
  outTree_->Branch("TnP_l2_mvaIdSpring15",TnP_l2_mvaIdSpring15,"TnP_l2_mvaIdSpring15[nTnP]/F");
  outTree_->Branch("TnP_l2_relIsoAfterFSR",TnP_l2_relIsoAfterFSR,"TnP_l2_relIsoAfterFSR[nTnP]/F");
  outTree_->Branch("TnP_l2_chargedHadIso03",TnP_l2_chargedHadIso03,"TnP_l2_chargedHadIso03[nTnP]/F");
  outTree_->Branch("TnP_l2_hasOwnFSR",TnP_l2_hasOwnFSR,"TnP_l2_hasOwnFSR[nTnP]/I");
  if(isMC()) {
    outTree_->Branch("TnP_l2_mcMatchId",TnP_l2_mcMatchId,"TnP_l2_mcMatchId[nTnP]/I");
    outTree_->Branch("TnP_l2_mcMatchAny",TnP_l2_mcMatchAny,"TnP_l2_mcMatchAny[nTnP]/I");
    outTree_->Branch("TnP_l2_mcPt",TnP_l2_mcPt,"TnP_l2_mcPt[nTnP]/F");
    outTree_->Branch("TnP_l2_mcPt1",TnP_l2_mcPt1,"TnP_l2_mcPt1[nTnP]/F");
  }
  outTree_->Branch("TnP_l2_hlt1L",TnP_l2_hlt1L,"TnP_l2_hlt1L[nTnP]/I");
  outTree_->Branch("TnP_l2_p4WithFSR_pt",TnP_l2_p4WithFSR_pt,"TnP_l2_p4WithFSR_pt[nTnP]/F");
  outTree_->Branch("TnP_l2_p4WithFSR_eta",TnP_l2_p4WithFSR_eta,"TnP_l2_p4WithFSR_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_p4WithFSR_phi",TnP_l2_p4WithFSR_phi,"TnP_l2_p4WithFSR_phi[nTnP]/F");
  outTree_->Branch("TnP_l2_p4WithFSR_mass",TnP_l2_p4WithFSR_mass,"TnP_l2_p4WithFSR_mass[nTnP]/F");
}


ZTnpAnalyzer::~ZTnpAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   

}


//
// member functions
//

// ------------ edm::Service<TFileService> outFile;
outTree_ = outFile->make<TTree>("zmmtree","muTnP tree for iso");method called for each event  ------------
void
ZTnpAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<double> fixedGridRhoFastjetAll;
   iEvent.getByToken(fixedGridRhoFastjetAllToken_,fixedGridRhoFastjetAll);
   const double evrho = *fixedGridRhoFastjetAll;

   edm::Handle<reco::VertexCollection> primaryVertices; 
   iEvent.getByToken(vertexToken_, primaryVertices);
   const reco::Vertex& vit = primaryVertices->front();

   const reco::Vertex *PV = 0;
   for (unsigned int i=0; i<primaryVertices->size(); i++) {
        PV = &(primaryVertices->at(i));        
        if (PV->chi2()==0 && PV->ndof()==0) continue;
        if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
        break;
   } 

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
  
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   

   PhysicsObjectSelector pObjsel(verbosity_);
   pObjsel.selectObjects(muons, electrons, pfs, jets, vit, evrho);

   const std::vector<std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>& tightMuon_ = *(pObjsel.tightSIPMuFSRpairVec());
   if(tightMuon_.size() < 2)         continue;
   const std::vector<pat::PackedCandidate>& fsrPVec_ = *(pObjsel.FSRVec());
   getZtnPpair(tightMuon_, fsrPVec_);
   outTree_->Fill();
   
}

void 
ZTnpAnalyze::getZtnPpair(const std::vector<std::pair< pat::Muon, std::vector<pat::PackedCandidate> >>& tightMuon,
                         const std::vector<pat::PackedCandidate>& fsrPVec) {
  int n = 0;
  std::vector<vhtm::ZtnP>  zvec;
  std::string hltPaths[] = {"HLT_IsoMu20_v", "HLT_IsoTkMu20_v", "HLT_IsoMu22_v", "HLT_IsoTkMu22_v", "HLT_IsoMu24_v", "HLT_IsoTkMu24_v"};
  for (unsigned int i = 0; i < tightMuon_.size(); ++i) {
    const auto& tagmu = tightMuon_[i].first;
    TLorentzVector tagP4 = HZZ4lUtil::getP4(tagmu);
    TLorentzVector tagP4wfsr = tagP4;
    bool taghasfsr = false;
    if(!tightMuon_[i].second.empty()) { 
      tagP4wfsr += HZZ4lUtil::getP4(tightMuon_[i].second.at(0));
      taghasfsr = true;
    }
    if(tagP4.Pt(edm::Service<TFileService> outFile;
outTree_ = outFile->make<TTree>("zmmtree","muTnP tree for iso");) <= 20.)    continue;  
    double tagiso = LeptonIsoCalculator::computeMuonReliso(tagmu, fsrPVec_);   
    if (tagiso >= 0.35) continue; // it is already scaled by pT
    bool tagMuhlt = false;
    for(unsigned int hli = 0; hli < 6; hli++ ) {
      const pat::TriggerObjectStandAloneCollection tagmuHLTMatches = tagmu.triggerObjectMatchesByPath( hltPaths[i] );
      if(!tagmuHLTMatches.empty()) {
        tagMuhlt = true;
        break;
      }
    }

    /*
    bool trigMatched = matchTriggerObject(tObj_, tagP4, "HLT_IsoMu20_v", 0, 30,matchedidx) < 0.02
    || matchTriggerObject(tObj_, tagP4, "HLT_IsoTkMu20_v", 0, 30,matchedidx) < 0.02
    || matchTriggerObject(tObj_, tagP4, "HLT_IsoMu22_v", 0, 30,matchedidx) < 0.02
    || matchTriggerObject(tObj_, tagP4, "HLT_IsoTkMu22_v", 0, 30,matchedidx) < 0.02;
    if (!trigMatched) continue; // it is already scaled by pT
    */
    for (unsigned int j = i+1; j < tightMuon_.size(); ++j) {
      const auto& probemu = tightMuon_[j].first;
      // opposite charge
      if ((tagmu.charge + probemu.charge) != 0) continue; 
      //std::cout << "Tag & Probe pair found!!" << std::endl;
      TLorentzVector probeP4 = HZZ4lUtil::getP4(probemu); 
      TLorentzVector probeP4wfsr = probeP4;
      bool probehasfsr = false;
      if(!tightMuon_[j].second.empty()) {
        probeP4wfsr += HZZ4lUtil::getP4(tightMuon_[j].second.at(0));
        probehasfsr = true;
      }
      double probeiso = LeptonIsoCalculator::computeMuonReliso(probemu, fsrPVec_, 0.01, 0.3); 
      //compute TnP pair P4
      TLorentzVector tnpP4 = tagP4 + probeP4;
      TLorentzVector tnpP4wfsr = tagP4wfsr + probeP4wfsr;
      vhtm::ZtnP ztemp;
      ztemp.TnP_pt = tnpP4.Pt();
      //Fill TnP pair property
      ztemp.TnP_pt = tnpP4.Pt();
      ztemp.TnP_eta = tnpP4.Eta();   
      ztemp.TnP_phi = tnpP4.Phi();   
      ztemp.TnP_mass = tnpP4wfsr.M();   //with fsr
      ztemp.TnP_hasFSR = taghasfsr || probehasfsr;   
      ztemp.TnP_mll = tnpP4.M();   //without fsr
      //tag property
      ztemp.TnP_l1_pdgId = (tagmu.charge > 0) ? -13 : 13;   
      ztemp.TnP_l1_pt = tagP4.Pt();   
      ztemp.TnP_l1_eta = tagP4.Eta();   
      ztemp.TnP_l1_phi = tagP4.Phi();   
      ztemp.TnP_l1_mass = tagP4.M();   
      ztemp.TnP_l1_charge = tagmu.charge;   
      ztemp.TnP_l1_tightId = 1;   
      ztemp.TnP_l1_looseId = 1;   
      ztemp.TnP_l1_dxy = tagmu.dxyPV;   
      ztemp.TnP_l1_dz  = tagmu.dzPV;   
      ztemp.TnP_l1_edxy = -1;   //not saved in our tree 
      ztemp.TnP_l1_edz = -1;//not saved in out tree
      ztemp.TnP_l1_ip3d = tagmu.dB3D;   //check
      ztemp.TnP_l1_sip3d = std::fabs(tagmu.dB3D/tagmu.edB3D);   //check
      ztemp.TnP_l1_ptErr = -1;//NS   
      ztemp.TnP_l1_lostHits = 0;//NS   
      ztemp.TnP_l1_trackerLayers = tagmu.trkHits;   
      ztemp.TnP_l1_pixelLayers = tagmu.pixHits;   
      ztemp.TnP_l1_etaSc = -1;   
      ztemp.TnP_l1_isGap = -1;   
      ztemp.TnP_l1_r9 = -1;   
      ztemp.TnP_l1_convVeto = -1;   
      ztemp.TnP_l1_mvaIdSpring15 = -1;;   
      ztemp.TnP_l1_relIsoAfterFSR = tagiso;   
      ztemp.TnP_l1_chargedHadIso03 = tagmu.pfChargedHadIsoR03;   
      ztemp.TnP_l1_hasOwnFSR  = taghasfsr; //doubt
      /*if(isMC()) {
        TLorentzVector genP4;
        int mid = 0;
        int genId = GenLevelMatching(tagP4, genObj_, genP4, mid);
        if(genId == -1) {
          //std::cout << "NOGEN" << std::endl;
          ztemp.TnP_l1_mcMatchId = 0;
          ztemp.TnP_l1_mcMatchAny = 0;
          ztemp.TnP_l1_mcPt = 0.;
          ztemp.TnP_l1_mcPt1 = 0.;
        } else {
          //if(abs(mid) != 23)   std::cout << "MIDMatched=" << mid << std::endl;
          ztemp.TnP_l1_mcMatchId = mid;
          ztemp.TnP_l1_mcMatchAny = 1;
          ztemp.TnP_l1_mcPt = genP4.Pt();
          ztemp.TnP_l1_mcPt1 = genP4.Pt();
        }
      }*/
      ztemp.TnP_l1_hlt1L = tagMuhlt; //doubt  t 
      ztemp.TnP_l1_p4WithFSR_pt = tagP4wfsr.Pt();   
      ztemp.TnP_l1_p4WithFSR_eta = tagP4wfsr.Eta();   
      ztemp.TnP_l1_p4WithFSR_phi = tagP4wfsr.Phi();   
      ztemp.TnP_l1_p4WithFSR_mass = tagP4wfsr.M();   
      //probe property
      ztemp.TnP_l2_pdgId = (probemu.charge > 0) ? -13 : 13;   
      ztemp.TnP_l2_pt = probeP4.Pt();   
      ztemp.TnP_l2_eta = probeP4.Eta();   
      ztemp.TnP_l2_phi = probeP4.Phi();   
      ztemp.TnP_l2_mass = probeP4.M();   
      ztemp.TnP_l2_charge = probemu.charge;   
      ztemp.TnP_l2_tightId = 1;   
      ztemp.TnP_l2_looseId = 1;   
      ztemp.TnP_l2_dxy = probemu.dxyPV;   
      ztemp.TnP_l2_dz  = probemu.dzPV;   
      ztemp.TnP_l2_edxy = -1;   //correct? 
      ztemp.TnP_l2_edz = -1;
      ztemp.TnP_l2_ip3d = tagmu.dB3D;   //check
      ztemp.TnP_l2_sip3d = std::fabs(tagmu.dB3D/tagmu.edB3D);   //check
      ztemp.TnP_l2_ptErr = -1;   
      ztemp.TnP_l2_lostHits = 0;   
      ztemp.TnP_l2_trackerLayers = probemu.trkHits;   
      ztemp.TnP_l2_pixelLayers = probemu.pixHits;   
      ztemp.TnP_l2_etaSc = -1;   
      ztemp.TnP_l2_isGap = -1;   
      ztemp.TnP_l2_r9 = -1;   
      ztemp.TnP_l2_convVeto = -1;   
      ztemp.TnP_l2_mvaIdSpring15 = -1;;   
      ztemp.TnP_l2_relIsoAfterFSR = probeiso;
      //if(probeiso > 0.35) std::cout << "FAILINGPROBE" << std::endl;      
      ztemp.TnP_l2_chargedHadIso03 = probemu.pfChargedHadIsoR03;   
      ztemp.TnP_l2_hasOwnFSR  = probehasfsr; //doubt
      /*if(isMC()) {
        TLorentzVector genP4;
        int mid = 0;
        int genId = GenLevelMatching(probeP4, genObj_, genP4, mid);
        if(genId == -1) {
          //std::cout << "NOGEN" << std::endl;
          ztemp.TnP_l2_mcMatchId = 0;
          ztemp.TnP_l2_mcMatchAny = 0;
          ztemp.TnP_l2_mcPt = 0.;
          ztemp.TnP_l2_mcPt1 = 0.;
        } else {
          //if(abs(mid) != 23)   std::cout << "MIDMatched=" << mid << std::endl;
          ztemp.TnP_l2_mcMatchId = mid;
          ztemp.TnP_l2_mcMatchAny = 1;
          ztemp.TnP_l2_mcPt = genP4.Pt();
          ztemp.TnP_l2_mcPt1 = genP4.Pt();
        }
      }*/
      /*ztemp.TnP_l2_hlt1L = matchTriggerObject(tObj_, probeP4, "HLT_IsoMu20_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoTkMu20_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoMu22_v", 0, 30,matchedidx) < 0.02
                           || matchTriggerObject(tObj_, probeP4, "HLT_IsoTkMu22_v", 0, 30,matchedidx) < 0.02;
      */ 
      ztemp.TnP_l2_p4WithFSR_pt = probeP4wfsr.Pt();   
      ztemp.TnP_l2_p4WithFSR_eta = probeP4wfsr.Eta();   
      ztemp.TnP_l2_p4WithFSR_phi = probeP4wfsr.Phi();   
      ztemp.TnP_l2_p4WithFSR_mass = probeP4wfsr.M(); 
      n++;
      zvec.push_back(ztemp);
    }
  }
  nTnP = n;
  for(int i = 0; i<n; i++) {
    if(i>=kMaxTnP)   continue;
    TnP_pt[i] = zvec[i].TnP_pt;
    
    TnP_eta[i] = zvec[i].TnP_eta;   
    TnP_phi[i] = zvec[i].TnP_phi;  
    TnP_mass[i] = zvec[i].TnP_mass;   
    TnP_hasFSR[i] = zvec[i].TnP_hasFSR;   
    TnP_mll[i] = zvec[i].TnP_mll;   
    TnP_l1_pdgId[i] = zvec[i].TnP_l1_pdgId;   
    TnP_l1_pt[i] = zvec[i].TnP_l1_pt;   
    TnP_l1_eta[i] = zvec[i].TnP_l1_eta;   
    TnP_l1_phi[i] = zvec[i].TnP_l1_phi;   
    TnP_l1_mass[i] = zvec[i].TnP_l1_mass;   
    TnP_l1_charge[i] = zvec[i].TnP_l1_charge;   
    TnP_l1_tightId[i] = zvec[i].TnP_l1_tightId;   
    TnP_l1_looseId[i] = zvec[i].TnP_l1_looseId;   
    TnP_l1_dxy[i] = zvec[i].TnP_l1_dxy;   
    TnP_l1_dz[i] = zvec[i].TnP_l1_dz;   
    TnP_l1_edxy[i] = zvec[i].TnP_l1_edxy;   
    TnP_l1_edz[i] = zvec[i].TnP_l1_edz;   
    TnP_l1_ip3d[i] = zvec[i].TnP_l1_ip3d;   
    TnP_l1_sip3d[i] = zvec[i].TnP_l1_sip3d;   
    TnP_l1_ptErr[i] = zvec[i].TnP_l1_ptErr;   
    TnP_l1_lostHits[i] = zvec[i].TnP_l1_lostHits;   
    TnP_l1_trackerLayers[i] = zvec[i].TnP_l1_trackerLayers;   
    TnP_l1_pixelLayers[i] = zvec[i].TnP_l1_pixelLayers;   
    TnP_l1_etaSc[i] = zvec[i].TnP_l1_etaSc;   
    TnP_l1_isGap[i] = zvec[i].TnP_l1_isGap;   
    TnP_l1_r9[i] = zvec[i].TnP_l1_r9;   
    TnP_l1_convVeto[i] = zvec[i].TnP_l1_convVeto;   
    TnP_l1_mvaIdSpring15[i] = zvec[i].TnP_l1_mvaIdSpring15;   
    TnP_l1_relIsoAfterFSR[i] = zvec[i].TnP_l1_relIsoAfterFSR;   
    TnP_l1_chargedHadIso03[i] = zvec[i].TnP_l1_chargedHadIso03;   
    TnP_l1_hasOwnFSR[i] = zvec[i].TnP_l1_hasOwnFSR;  
    if(isMC()) {
      TnP_l1_mcMatchId[i] = zvec[i].TnP_l1_mcMatchId;
      TnP_l1_mcMatchAny[i] = zvec[i].TnP_l1_mcMatchAny;
      TnP_l1_mcPt[i] = zvec[i].TnP_l1_mcPt;
      TnP_l1_mcPt1[i] = zvec[i].TnP_l1_mcPt1;
    } 
    TnP_l1_hlt1L[i] = zvec[i].TnP_l1_hlt1L;   
    TnP_l1_p4WithFSR_pt[i] = zvec[i].TnP_l1_p4WithFSR_pt;   
    TnP_l1_p4WithFSR_eta[i] = zvec[i].TnP_l1_p4WithFSR_eta;   
    TnP_l1_p4WithFSR_phi[i] = zvec[i].TnP_l1_p4WithFSR_phi;   
    TnP_l1_p4WithFSR_mass[i] = zvec[i].TnP_l1_p4WithFSR_mass;   
    TnP_l2_pdgId[i] = zvec[i].TnP_l2_pdgId;   
    TnP_l2_pt[i] = zvec[i].TnP_l2_pt;   
    TnP_l2_eta[i] = zvec[i].TnP_l2_eta;   
    TnP_l2_phi[i] = zvec[i].TnP_l2_phi;   
    TnP_l2_mass[i] = zvec[i].TnP_l2_mass;   
    TnP_l2_charge[i] = zvec[i].TnP_l2_charge;   
    TnP_l2_tightId[i] = zvec[i].TnP_l2_tightId;   
    TnP_l2_looseId[i] = zvec[i].TnP_l2_looseId;   
    TnP_l2_dxy[i] = zvec[i].TnP_l2_dxy;   
    TnP_l2_dz[i] = zvec[i].TnP_l2_dz;   
    TnP_l2_edxy[i] = zvec[i].TnP_l2_edxy;   
    TnP_l2_edz[i] = zvec[i].TnP_l2_edz;   
    TnP_l2_ip3d[i] = zvec[i].TnP_l2_ip3d;   
    TnP_l2_sip3d[i] = zvec[i].TnP_l2_sip3d;   
    TnP_l2_ptErr[i] = zvec[i].TnP_l2_ptErr;   
    TnP_l2_lostHits[i] = zvec[i].TnP_l2_lostHits;   
    TnP_l2_trackerLayers[i] = zvec[i].TnP_l2_trackerLayers;   
    TnP_l2_pixelLayers[i] = zvec[i].TnP_l2_pixelLayers;   
    TnP_l2_etaSc[i] = zvec[i].TnP_l2_etaSc;   
    TnP_l2_isGap[i] = zvec[i].TnP_l2_isGap;   
    TnP_l2_r9[i] = zvec[i].TnP_l2_r9;   
    TnP_l2_convVeto[i] = zvec[i].TnP_l2_convVeto;   
    TnP_l2_mvaIdSpring15[i] = zvec[i].TnP_l2_mvaIdSpring15;   
    TnP_l2_relIsoAfterFSR[i] = zvec[i].TnP_l2_relIsoAfterFSR;   
    TnP_l2_chargedHadIso03[i] = zvec[i].TnP_l2_chargedHadIso03;   
    TnP_l2_hasOwnFSR[i] = zvec[i].TnP_l2_hasOwnFSR;   
    if(isMC()) {
      TnP_l2_mcMatchId[i] = zvec[i].TnP_l2_mcMatchId;
      TnP_l2_mcMatchAny[i] = zvec[i].TnP_l2_mcMatchAny;
      TnP_l2_mcPt[i] = zvec[i].TnP_l2_mcPt;
      TnP_l2_mcPt1[i] = zvec[i].TnP_l2_mcPt1; 
    }
    TnP_l2_hlt1L[i] = zvec[i].TnP_l2_hlt1L;   
    TnP_l2_p4WithFSR_pt[i] = zvec[i].TnP_l2_p4WithFSR_pt;   
    TnP_l2_p4WithFSR_eta[i] = zvec[i].TnP_l2_p4WithFSR_eta;   
    TnP_l2_p4WithFSR_phi[i] = zvec[i].TnP_l2_p4WithFSR_phi;   
    TnP_l2_p4WithFSR_mass[i] = zvec[i].TnP_l2_p4WithFSR_mass;
  }
  //std::cout << "NTnP=" << n << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZTnpAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZTnpAnalyzer::endJob() 
{
  std::cout << "Number of events with atleast 4 tight leptons=" << nEvtswith4Lep_ << std::endl;
  std::cout << "Number of events with atleast 2 Zcandidates=" << nEvtwith2Z_ << std::endl;
  std::cout << "Number of events with atleast 1 ZZcandidate=" << nEvtwithZZ_ << std::endl;
  std::cout << "Number of events with atleast 1 4mu candidate=" << nEvtwith4m_ << std::endl;
  std::cout << "Number of events with atleast 1 4e candidate=" << nEvtwith4e_ << std::endl;
  std::cout << "Number of events with atleast 1 2e2m candidate=" << nEvtwith2e2m_ << std::endl;
  outTree_->AutoSave();
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
