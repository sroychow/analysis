#include "DataFormats/Common/interface/Wrapper.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

#include <vector>
#include<map>
#include<string>
namespace {
  struct dictionary {
    vhtm::Event rv1;
    edm::Wrapper<vhtm::Event> wrv1;
    vhtm::GenEvent rv2;
    edm::Wrapper<vhtm::GenEvent> wrv2;
    vhtm::Electron rv3;
    edm::Wrapper<vhtm::Electron> wrv3;
    vhtm::GenParticle rv4;
    edm::Wrapper<vhtm::GenParticle> wrv4;
    vhtm::GenJet rv5;
    edm::Wrapper<vhtm::GenJet> wrv5;
    vhtm::GenMET rv6;
    edm::Wrapper<vhtm::GenMET> wrv6;
    vhtm::MET rv7;
    edm::Wrapper<vhtm::MET> wrv7;
    vhtm::Tau rv8;
    edm::Wrapper<vhtm::Tau> wrv8;
    vhtm::Muon rv9;
    edm::Wrapper<vhtm::Muon> wrv9;
    vhtm::Jet rva;
    edm::Wrapper<vhtm::Jet> wrva;
    vhtm::Vertex rvb;
    edm::Wrapper<vhtm::Vertex> wrvb;
    vhtm::TriggerObject rvd;
    edm::Wrapper<vhtm::TriggerObject> wrvd;
    vhtm::Candidate rve;
    edm::Wrapper<vhtm::Candidate> wrve;
    vhtm::Photon rvf;
    edm::Wrapper<vhtm::Photon> wrvf;
    vhtm::PackedPFCandidate rvg;
    edm::Wrapper<vhtm::PackedPFCandidate> wrvg;

    std::vector<vhtm::Event> vrvc;
    std::vector<vhtm::GenEvent> vrvd;
    std::vector<vhtm::Electron> vrv1;
    std::vector<vhtm::GenParticle> vrv2;
    std::vector<vhtm::GenJet> vrv3;
    std::vector<vhtm::GenMET> vrv4;
    std::vector<vhtm::MET> vrvb;    
    std::vector<vhtm::Tau> vrv5;
    std::vector<vhtm::Muon> vrv6;
    std::vector<vhtm::Jet> vrv7;
    std::vector<vhtm::Vertex> vrv8;
    std::vector<vhtm::TriggerObject> vrva;
    std::vector<vhtm::Candidate> vrve;
    std::vector<vhtm::Photon> vrvf;
    std::vector<vhtm::PackedPFCandidate> vrvg;
    
    std::vector<double> vrvl;
    std::map<double, std::vector<double> > vrvm;
    std::map< std::string, std::vector<double> > vrvn;
  };
}
