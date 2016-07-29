import FWCore.ParameterSet.Config as cms

bbtoDijet = cms.EDAnalyzer('bbtoDijetAnalyzer',
    processName     = cms.string("TEST"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD", "", "TEST"),
    triggerResults  = cms.InputTag("TriggerResults::TEST"),
    caloJets        = cms.InputTag("ak4CaloJets::RECO"),
    pfJets          = cms.InputTag("ak4PFJets::RECO"),
    bTagCSVOnline   = cms.InputTag("hltCombinedSecondaryVertexBJetTagsCalo"),
    bTagCSVOffline  = cms.InputTag("pfCombinedSecondaryVertexV2BJetTags")
)
