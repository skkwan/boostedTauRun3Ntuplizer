import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TEventDisplayGenerator",
                                  rctSource               = cms.InputTag("simRctDigis"),
                                  #gctTauJetsSource        = cms.InputTag("simCaloStage1LegacyFormatDigis","tauJets"),
                                  #gctIsoTauJetsSource     = cms.InputTag("simCaloStage1LegacyFormatDigis","isoTauJets"),
                                  #l1ExtraJetSource        = cms.InputTag("uct2016EmulatorDigis","Central"),
                                  #l1ExtraTauSource        = cms.InputTag("uct2016EmulatorDigis","Tau"),
                                  #l1ExtraIsoTauSource     = cms.InputTag("uct2016EmulatorDigis","IsoTau"),
                                  #recoTaus                 = cms.InputTag("hpsPFTauProducer"),
                                  #slimmedTaus              = cms.InputTag("slimmedTaus"),
                                  #recoJets                 = cms.InputTag("slimmedJets"),
                                  #remove all possible muons
                                  #recoTauDiscriminatorIso = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
                                  #recoTauDiscriminatorMu  = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
                                  vertices                = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  folderName              = cms.untracked.string("firstFolder"),
                                  recoPtCut               = cms.double(4),
                                  #packedPfCands           = cms.InputTag("packedPFCandidates"),
                                  #pfCands                 = cms.InputTag("particleFlow"),
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis'),
                                  ecalCaloClusters        = cms.InputTag('particleFlowSuperClusterECAL','particleFlowBasicClusterECALBarrel','RECO')
)
