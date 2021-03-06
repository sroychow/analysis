import FWCore.ParameterSet.Config as cms

triggerBlock = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.untracked.int32(0),
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
  ),
  l1Tag = cms.InputTag('gtDigis'),
  hltTag = cms.InputTag('TriggerResults','','HLT')
)

triggerBlockZTnP = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.untracked.int32(0),
  hltPathsOfInterest = cms.vstring ( 
 	'HLT_IsoMu20_v*',
  	'HLT_IsoTkMu20_v*',
        'HLT_IsoMu22_v*',
        'HLT_IsoTkMu22_v*'
  ),
  l1Tag = cms.InputTag('gtDigis'),
  hltTag = cms.InputTag('TriggerResults','','HLT')
)
