import FWCore.ParameterSet.Config as cms

from AnalysisSpace.TreeMaker.EventBlock_cfi import eventBlock
from AnalysisSpace.TreeMaker.VertexBlock_cfi import vertexBlock
from AnalysisSpace.TreeMaker.ElectronBlock_bx25_cfi import electronBlock
from AnalysisSpace.TreeMaker.CalibratedElectronBlock_bx25_cfi import calibelectronBlock
from AnalysisSpace.TreeMaker.MuonBlock_cfi import muonBlock
from AnalysisSpace.TreeMaker.TauBlock_cfi import tauBlock
from AnalysisSpace.TreeMaker.JetBlock_cfi import jetBlock
from AnalysisSpace.TreeMaker.METBlock_cfi import metBlock
from AnalysisSpace.TreeMaker.TriggerBlock_cfi import triggerBlock
from AnalysisSpace.TreeMaker.TriggerObjectBlock_cfi import triggerObjectBlock
from AnalysisSpace.TreeMaker.PackedPFCandidateBlock_cfi import packedPFCandidateBlock

treeContentSequence = cms.Sequence(
   eventBlock
 + vertexBlock
 + electronBlock
 + calibelectronBlock
 + muonBlock
# + tauBlock
 + jetBlock
 + metBlock
 + triggerBlock
 + triggerObjectBlock
 + packedPFCandidateBlock
)
