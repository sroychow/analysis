#include <iostream>
#include <algorithm>
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "AnalysisSpace/OnlineAnalysis/plugins/FourLeptonSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/PhysicsObjectSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"
#include "AnalysisSpace/OnlineAnalysis/src/ZCandidateSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/ZZCandidateSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/KinZrefitter.h"
// constructors and destructor
//
FourLeptonSelector::FourLeptonSelector(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  syncFilename_(iConfig.getUntrackedParameter<std::string>("syncFilename")),
  fixedGridRhoFastjetAllToken_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllTag",edm::InputTag("fixedGridRhoFastjetAll")))),  
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcandSrc"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   nEvtswith4Lep_ = 0;
   nEvtwith2Z_ = 0;  
   nEvtwithZZ_ = 0;
   nEvtwith2e2m_ = 0;
   nEvtwith4m_ = 0;
   nEvtwith4e_ = 0; 
   //_(consumes<>(iConfig.getParameter<edm::InputTag>())),
   mela_ = new Mela(13.0, 125.0, TVar::SILENT);
   mela_->setCandidateDecayMode(TVar::CandidateDecay_ZZ);  
   syncDumpf_.open(syncFilename_.c_str(), ios::out);
   if (!syncDumpf_) std::cerr << "Output File: " << syncFilename_ << " could not be opened!" << endl;
}


FourLeptonSelector::~FourLeptonSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FourLeptonSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   edm::Handle<double> fixedGridRhoFastjetAll;
   iEvent.getByToken(fixedGridRhoFastjetAllToken_,fixedGridRhoFastjetAll);
   const double evrho = *fixedGridRhoFastjetAll;

   edm::Handle<reco::VertexCollection> primaryVertices; 
   iEvent.getByToken(vertexToken_, primaryVertices);
   const reco::Vertex& vit = primaryVertices->front();

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
  
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   
   PhysicsObjectSelector pObjsel(verbosity_);
   pObjsel.selectObjects(muons, electrons, pfs, vit, evrho);

   std::vector<vhtm::IsoElectron>*  iele = pObjsel.tightIsoEleFSRpairVec();
   std::vector<vhtm::IsoMuon>* imu = pObjsel.tightIsoMuFSRpairVec();
   //ignore event if tight electron + muon < 4
   if(iele->size() + imu->size() >= 4)  nEvtswith4Lep_++;
   ZCandidateSelector zSel;
   zSel.selectZcandidates(*imu, *iele);
   std::vector<vhtm::Zee>* zeeV = zSel.getZeeVec();
   std::vector<vhtm::Zmumu>* zmmV = zSel.getZmumuVec();  
   if(zeeV->size() + zmmV->size() >= 2)   nEvtwith2Z_++;
   ZZCandidateSelector zzSel;
   zzSel.selectZZcandidates(*zeeV, *zmmV);
   std::vector<vhtm::ZZcandidate> zzcandV = (*zzSel.getZZVec());
   bool m4 = false;
   bool e4 = false;
   bool m2e2 = false;
   bool foundbestZ = false;
 
  if(zzcandV.size() >= 1)  foundbestZ = true;

   vhtm::ZZcandidate bestZZcand;
   if(zzcandV.size() == 1)   {
     nEvtwithZZ_++;
     bestZZcand = zzcandV[0];
     KDCalculator kdc;
     kdc.computeKDforbestZZ(mela_, bestZZcand);
   } 
   else if(zzcandV.size() > 1)   {
     nEvtwithZZ_++;
     KDCalculator kdc;
     for(auto& zz : zzcandV){
       kdc.computeKDforbestZZ(mela_, zz);
     }
     //If more than one ZZ candidate survives, choose the one with the highest value of P_sig/P_bkg (p0plus_VAJHU/bkg_VAMCFM). 
     //To choose among candidates that use the same 4 leptons (alternate pairings in 4mu/4e), 
     //always choose the one with mZ1 closest to nominal mZ. 
     std::sort(zzcandV.begin(), zzcandV.end(), HZZ4lUtil::ZZkdComparator());
     bestZZcand = zzcandV[0];
     for(unsigned int i = 1; i < zzcandV.size(); i++) {
       if(bestZZcand.Dbkgkin() > zzcandV[i].Dbkgkin()) {
         //foundbestZ = true;//no need to check others   
         break;
       }
       else if(bestZZcand.Dbkgkin() == zzcandV[i].Dbkgkin()){
         if(bestZZcand.mZ1 > zzcandV[i].mZ1) continue;          
         else if(bestZZcand.mZ1 == zzcandV[i].mZ1) {//for same Z1, choose the one with Z2 higset sumptlepton
           if(bestZZcand.Z2lepPtavg > zzcandV[i].Z2lepPtavg)  continue;
           else bestZZcand = zzcandV[i];
         } else bestZZcand = zzcandV[i];
       }  
     }
   }
   if(foundbestZ) {
     if(bestZZcand.flavour == HZZ4lUtil::ZZType::mmmm) m4 = true;
     else if(bestZZcand.flavour == HZZ4lUtil::ZZType::eeee) e4 = true;
     else if(bestZZcand.flavour == HZZ4lUtil::ZZType::eemm || bestZZcand.flavour == HZZ4lUtil::ZZType::mmee)  m2e2 = true;
     if(m2e2) nEvtwith2e2m_++;
     if(m4) nEvtwith4m_++;
     if(e4) nEvtwith4e_++;
     KinZrefitter kinzr(isMC_);
     kinzr.doZmassrefit(bestZZcand);
     HZZ4lUtil::syncDumper(iEvent.id().run(), iEvent.id().luminosityBlock(),iEvent.id().event(), bestZZcand, syncDumpf_);
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
FourLeptonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FourLeptonSelector::endJob() 
{
  std::cout << "Number of events with atleast 4 tight leptons=" << nEvtswith4Lep_ << std::endl;
  std::cout << "Number of events with atleast 2 Zcandidates=" << nEvtwith2Z_ << std::endl;
  std::cout << "Number of events with atleast 1 ZZcandidate=" << nEvtwithZZ_ << std::endl;
  std::cout << "Number of events with atleast 1 4mu candidate=" << nEvtwith4m_ << std::endl;
  std::cout << "Number of events with atleast 1 4e candidate=" << nEvtwith4e_ << std::endl;
  std::cout << "Number of events with atleast 1 2e2m candidate=" << nEvtwith2e2m_ << std::endl;
  syncDumpf_.close();
  delete mela_;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FourLeptonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FourLeptonSelector);
