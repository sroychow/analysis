#include <iostream>

#include "TTree.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "Utilities/General/interface/FileInPath.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/TauBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

namespace tb
{
  template<typename T>
  bool isValidRef(const edm::Ref<T>& ref) {
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
  }
}
TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  tauTag_(iConfig.getUntrackedParameter<edm::InputTag>("patTauSrc", edm::InputTag("selectedPatTaus"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", true)),
  tauToken_(consumes<pat::TauCollection>(tauTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_))
{
  produces<std::vector<vhtm::Tau>>("vhtmTauVector").setBranchAlias("vhtmTauVector");
}
TauBlock::~TauBlock() {}
void TauBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Tau>();
  tree->Branch("Tau", "std::vector<vhtm::Tau>", &list_, 32000, -1);
  tree->Branch("nTau", &fnTau_, "fnTau_/I");
}
void TauBlock::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnTau_ = 0;

  edm::Handle<pat::TauCollection> taus;
  bool found = iEvent.getByToken(tauToken_, taus);
  
  if (found && taus.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    if (bsCorr_)
      iEvent.getByToken(bsToken_, beamSpot);

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    for (const pat::Tau& v: *taus) {
      if (list_->size() == kMaxTau_) {
        edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << list_->size();
        break;
      }
      vhtm::Tau tau;

      // Store Tau variables
      tau.eta    = v.eta();
      tau.phi    = v.phi();
      tau.pt     = v.pt();
      tau.energy = v.energy();
      tau.charge = v.charge();
#if 0
      if (v.leadChargedHadrCand().isNonnull() && tb::isValidRef(v.leadChargedHadrCand()->trackRef())) {
        reco::TrackRef trk = v.leadChargedHadrCand()->trackRef();
        tau.leadTrkPt      = trk->pt();

        if (primaryVertices.isValid()) {
          edm::LogInfo("TauBlock") << "Total # Primary Vertices: " 
                                   << primaryVertices->size();

          // IP of leadChargedHadrCand wrt event PV
          auto vit = primaryVertices->begin(); // Highest sumPt vertex
          tau.dxyPV = trk->dxy(vit->position());
          tau.dzPV  = trk->dz(vit->position());

          // IP of leadChargedHadrCand wrt closest PV
          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDz = 9999.;
          double vertexDxy = 9999.;
          for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
            double dxy = trk->dxy(vit->position());
            double dz  = trk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDxy = dxy;
              vertexDz = dz;
            }
          }
          tau.vtxIndex = indexVtx;
          tau.vtxDxy   = vertexDxy;
          tau.vtxDz    = vertexDz;
        }
        else {
          edm::LogError("TauBlock") << "Error >> Failed to get VertexCollection for label: "
                                    << vertexTag_;
        }
      }
#endif
      // Leading particle pT
      tau.leadChargedParticlePt = 
           v.leadChargedHadrCand().isNonnull() ? v.leadChargedHadrCand()->pt(): 0.;
      tau.leadNeutralParticlePt = 
           v.leadNeutralCand().isNonnull() ? v.leadNeutralCand()->et(): 0.;
      tau.leadParticlePt = 
           v.leadCand().isNonnull() ? v.leadCand()->et(): 0.;

      // Signal Constituents
      // Charged hadrons
      for (const reco::CandidatePtr& iCand: v.signalChargedHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.sigChHadList.push_back(c);
      }
      // Neutral hadrons
      for (const reco::CandidatePtr& iCand: v.signalNeutrHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.sigNeHadList.push_back(c);
      }
      // Photons
      for (const reco::CandidatePtr& iCand: v.signalGammaCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.sigGammaList.push_back(c);
      }
      // Isolation Constituents
      // Charged hadrons
      double isoPtSumChHad = 0;
      for (const reco::CandidatePtr& iCand: v.isolationChargedHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.isoChHadList.push_back(c);
        isoPtSumChHad += cand.pt();
      }
      // Neutral hadrons
      double isoPtSumNeHad = 0;
      for (const reco::CandidatePtr& iCand: v.isolationNeutrHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.isoNeHadList.push_back(c);
        isoPtSumNeHad += cand.pt();
      }
      // Photons
      double isoPtSumGamma = 0;
      for (const reco::CandidatePtr& iCand: v.isolationGammaCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.isoGammaList.push_back(c);
        isoPtSumGamma += cand.pt();
      }
      // PtSum
      tau.ptSumChargedHadronsIsoCone = isoPtSumChHad;
      tau.ptSumNeutralHadronsIsoCone = isoPtSumNeHad;
      tau.ptSumPhotonsIsoCone        = isoPtSumGamma;

      // tau id. discriminators
      tau.decayModeFinding = v.tauID("decayModeFinding");
      tau.decayModeFindingNewDMs = v.tauID("decayModeFindingNewDMs");
//      tau.decayModeFindingOldDMs = v.tauID("decayModeFindingOldDMs");

      // discriminators against muons
      //tau.againstMuonLoose  = v.tauID("againstMuonLoose");
      //tau.againstMuonMedium = v.tauID("againstMuonMedium");
      //tau.againstMuonTight  = v.tauID("againstMuonTight");

      tau.againstMuonLoose3  = v.tauID("againstMuonLoose3");
      tau.againstMuonTight3 = v.tauID("againstMuonTight3");

      // discriminators against electrons
      //tau.againstElectronLoose  = v.tauID("againstElectronLoose");
      //tau.againstElectronMedium = v.tauID("againstElectronMedium");
      //tau.againstElectronTight  = v.tauID("againstElectronTight");

      tau.againstElectronLooseMVA5  = v.tauID("againstElectronLooseMVA5");
      tau.againstElectronMediumMVA5 = v.tauID("againstElectronMediumMVA5");
      tau.againstElectronTightMVA5  = v.tauID("againstElectronTightMVA5");

      // DB Corrected Isolation
      tau.byLooseCombinedIsolationDeltaBetaCorr3Hits  = v.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      tau.byMediumCombinedIsolationDeltaBetaCorr3Hits = v.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tau.byTightCombinedIsolationDeltaBetaCorr3Hits  = v.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      tau.byCombinedIsolationDeltaBetaCorrRaw3Hits    = v.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

      // Isolation variables suggested by Christian
      tau.chargedIsoPtSum = v.tauID("chargedIsoPtSum");
      tau.neutralIsoPtSum = v.tauID("neutralIsoPtSum");
      tau.puCorrPtSum     = v.tauID("puCorrPtSum");
#if 0
      std::cout << ">>> tauID(label) = value" << std::endl;
      for (const pat::Tau::IdPair& pa: v.tauIDs())
	std::cout << pa.first << "=" << pa.second << std::endl;
#endif
      // kinematic variables for PFJet associated to PFTau
#if 0
      if (tb::isValidRef(v.pfJetRef())) {
        const reco::PFJetRef &jtr = v.pfJetRef();
        tau.jetPt  = jtr->pt();
        tau.jetEta = jtr->eta();
        tau.jetPhi = jtr->phi();
      }
      // NEW quantities
      tau.emFraction = v.emFraction();
#endif

      // Vertex information
      const reco::Candidate::Point& vertex = v.vertex();
      tau.vx = vertex.x();
      tau.vy = vertex.y();
      tau.vz = vertex.z();

      tau.zvertex = v.vz(); // distance from the primary vertex
      tau.mass    = v.p4().M();
#if 0
      tau.ltsipt  = TMath::Abs(v.leadChargedHadrCandsignedSipt());
#endif
      // pat::Tau now has reference to the primary vertex also

      // add particle to the list
      list_->push_back(tau);
    }
    fnTau_ = list_->size();
    //put the vhtm collections in edm
    std::auto_ptr<std::vector<vhtm::Tau>> pv1(new std::vector<vhtm::Tau>(*list_));
    iEvent.put(pv1,"vhtmTauVector");
  }
  else {
    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: "
                              << tauTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);

