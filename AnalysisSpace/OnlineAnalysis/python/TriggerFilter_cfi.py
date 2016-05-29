import FWCore.ParameterSet.Config as cms

triggerFilter = cms.EDFilter("TriggerFilter",
  verbosity = cms.untracked.int32(0),
  l1InputTag = cms.untracked.InputTag('gtDigis'),
  hltInputTag = cms.untracked.InputTag('TriggerResults','','HLT'),
  hltPathsOfInterest = cms.vstring ( 
    'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*',
    'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
    'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',
    'HLT_TripleMu_12_10_5_v*',
    'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*',
    'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*',
    'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*',
    'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*',
    'HLT_Ele23_WPLoose_Gsf_v*' 
  )
)
