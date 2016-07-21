import FWCore.ParameterSet.Config as cms

hltana = cms.EDAnalyzer('SimpleHLTAnalyzer',
    hltProcess = cms.string('TEST'),
    # This is required to access basic trigger information -- i.e, which event were accepted by which triggers
    trigSummary = cms.InputTag('hltTriggerSummaryAOD', '', 'TEST'),
    trigResults = cms.InputTag('TriggerResults::TEST'),
    # For accessing offline trigger objects
    pfJets = cms.InputTag('ak4PFJets::RECO'),
    bTagsCSVOffline = cms.InputTag('pfCombinedSecondaryVertexV2BJetTags'),
    bTagsCSVOnline  = cms.InputTag('hltCombinedSecondaryVertexBJetTagsCalo')
)
