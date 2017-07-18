#include <iostream>
#include <algorithm>
#include "AnalysisSpace/OnlineAnalysis/interface/OnlineObjects.h"
#include "AnalysisSpace/OnlineAnalysis/plugins/SelectedLeptonProducer.h"
#include "AnalysisSpace/OnlineAnalysis/src/PhysicsObjectSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/HZZ4lUtil.h"
#include "AnalysisSpace/OnlineAnalysis/src/ZCandidateSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/ZZCandidateSelector.h"
#include "AnalysisSpace/OnlineAnalysis/src/KinZrefitter.h"
// constructors and destructor
//
SelectedLeptonProducer::SelectedLeptonProducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  fixedGridRhoFastjetAllToken_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllTag",edm::InputTag("fixedGridRhoFastjetAll")))),  
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcandSrc"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc")))
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Electron> >();
   produces<std::vector<pat::Muon> >();
   produces<std::vector<pat::PackedCandidate> >();
   produces<std::vector<pat::JetCollection> >();
}


SelectedLeptonProducer::~SelectedLeptonProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SelectedLeptonProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
/*   std::cout << "PV" 
   << std::setprecision(2)
     << std::setw(8) << PV->ndof()
     << std::setw(8) << PV->position().Z() 
     << std::setw(8) << PV->chi2()
     << std::setw(8) << PV->position().Rho()
      << std::endl;
*/
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
  
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);

   PhysicsObjectSelector pObjsel(verbosity_);
   pObjsel.selectObjects(muons, electrons, pfs, jets, vit, evrho);

   std::vector<vhtm::IsoElectron>*  iele = pObjsel.tightIsoEleFSRpairVec();
   std::vector<vhtm::IsoMuon>* imu = pObjsel.tightIsoMuFSRpairVec();

}


// ------------ method called once each job just before starting event loop  ------------
void 
SelectedLeptonProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SelectedLeptonProducer::endJob() 
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
SelectedLeptonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SelectedLeptonProducer);
