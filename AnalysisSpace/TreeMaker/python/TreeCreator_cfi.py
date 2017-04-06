import FWCore.ParameterSet.Config as cms

treeCreator = cms.EDAnalyzer("TreeMakerModule",
  verbosity = cms.untracked.int32(1),
  createTree = cms.bool(True),
  treeName = cms.string("vhtree")
)
#online analysis tree
treeCreatorOnline = cms.EDAnalyzer("TreeMakerModule",
  verbosity = cms.untracked.int32(1),
  createTree = cms.bool(True),
  treeName = cms.string("onlineTree")
)
