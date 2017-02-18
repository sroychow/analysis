#include <iostream>
#include <algorithm>
#include <TRandom3.h>
#include "AnalysisSpace/OnlineAnalysis/plugins/JetQGproducer.h"
// JEC related
#include "PhysicsTools/PatAlgos/plugins/PATJetUpdater.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
 
//JER related
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


//
// constructors and destructor
//
JetQGproducer::JetQGproducer(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  ptCut_(iConfig.getParameter<double>("jetPtcut")),
  etaCut_(iConfig.getParameter<double>("jetEtacut")),
  muRhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("muRho"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  qgToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood")))
{
  produces<std::vector<pat::Jet> >(); 
}


JetQGproducer::~JetQGproducer()
{
}


// ------------ method called to produce the data  ------------
void
JetQGproducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Rho Correction
  edm::Handle<double> evrmu;
  iEvent.getByToken(muRhoToken_, evrmu);
  double muRho = *evrmu;

  // JER
  //see example https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#CMSSW_7_6_X_and_CMSSW_8_0X
  JME::JetResolution resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  JME::JetResolution resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");
  JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
  
  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(jetToken_, jets);

  //new collection to write
  std::vector<pat::Jet>*   qgJets = new std::vector<pat::Jet>();

  //check example on twiki
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood#CMSSW_8_0_20_X_recomendation
  edm::Handle<edm::ValueMap<float> > qgHandle;
  iEvent.getByToken(qgToken_, qgHandle);

  //Steps
  //1. Skim the jet collection. Pt > 30. , |eta| < 4.7
  //2. Add the q/g value from the valuemap
  //If MC: apply smearing
  for(auto jet = jets->begin();  jet != jets->end(); ++jet) {
    if(jet->pt() <= ptCut_)     continue;
    if(std::fabs(jet->eta()) >= etaCut_)    continue;
    edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jets, jet - jets->begin()));
    float qgLikelihood = (*qgHandle)[jetRef]; 
    //JER part;From HZZ4l twiki
    JME::JetParameters parameters;
    parameters.setJetPt(jet->pt());
    parameters.setJetEta(jet->eta());
    parameters.setRho(muRho);
    float relpterr = resolution_pt.getResolution(parameters);
    //float phierr = resolution_phi.getResolution(parameters);
    
    double jercorr = 1.0; double jercorrup = 1.0; double jercorrdn = 1.0;
    if (isMC_) {
      JME::JetParameters sf_parameters = {{JME::Binning::JetEta, jet->eta()}, {JME::Binning::Rho, muRho}};
      float factor = resolution_sf.getScaleFactor(sf_parameters);
      float factorup = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
      float factordn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);
      
      double pt_jer, pt_jerup, pt_jerdn;
      const reco::GenJet * genJet = jet->genJet();
      if (genJet && deltaR(jet->eta(),jet->phi(),genJet->eta(),genJet->phi())<0.2
	  && (std::abs(jet->pt()-genJet->pt())<3*relpterr*jet->pt()) ) {
	double gen_pt = genJet->pt();
	pt_jer = std::max(0.0,gen_pt+factor*(jet->pt()-gen_pt));
	pt_jerup = std::max(0.0,gen_pt+factorup*(jet->pt()-gen_pt));
	pt_jerdn = std::max(0.0,gen_pt+factordn*(jet->pt()-gen_pt));
      } else {
	TRandom3 rand;
	rand.SetSeed(std::abs(static_cast<int>(sin(jet->phi())*100000)));
	float smear = rand.Gaus(0,1.0);
	float sigma = sqrt(factor*factor-1.0)*relpterr*jet->pt();
	float sigmaup = sqrt(factorup*factorup-1.0)*relpterr*jet->pt();
	float sigmadn = sqrt(factordn*factordn-1.0)*relpterr*jet->pt();
	pt_jer = std::max(0.0,smear*sigma+jet->pt());
	pt_jerup = std::max(0.0,smear*sigmaup+jet->pt());
	pt_jerdn = std::max(0.0,smear*sigmadn+jet->pt());
      }
      
      jercorr = pt_jer/jet->pt();
      jercorrup = pt_jerup/jet->pt();
      jercorrdn = pt_jerdn/jet->pt();
    }
    
    TLorentzVector *jet_jer = new TLorentzVector(jercorr*jet->px(),jercorr*jet->py(),jercorr*jet->pz(),jercorr*jet->energy());
    TLorentzVector *jet_jerup = new TLorentzVector(jercorrup*jet->px(),jercorrup*jet->py(),jercorrup*jet->pz(),jercorrup*jet->energy());
    TLorentzVector *jet_jerdn = new TLorentzVector(jercorrdn*jet->px(),jercorrdn*jet->py(),jercorrdn*jet->pz(),jercorrdn*jet->energy());
    pat::Jet qJet(*jet);
    qJet.addUserFloat("qgLikelihood",qgLikelihood); 
    qJet.addUserData("jet_jerP4",*jet_jer);   
    qJet.addUserData("jet_jerupP4",*jet_jerup);   
    qJet.addUserData("jet_jerdnP4",*jet_jerdn);   
    qgJets->push_back(qJet);
  }
  
  std::auto_ptr<std::vector<pat::Jet> > ptr(qgJets);
  iEvent.put(ptr);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
JetQGproducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
JetQGproducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  JetQGproducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  JetQGproducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
JetQGproducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

 
// ------------ method called when ending the processing of a luminosity block  ------------

void
JetQGproducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetQGproducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetQGproducer);
