import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("Run3Ntuplizer",
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  recoJets                = cms.InputTag("slimmedJets"),
                                  recoJetsAK8             = cms.InputTag("slimmedJetsAK8"),
                                  miniTaus                = cms.InputTag("slimmedTaus"),
                                  genParticles     = cms.InputTag("genParticles", "", "HLT"),
                                  recoPtCut               = cms.double(10),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis'),
                                  l1UCTCentralJets        = cms.InputTag("uct2016EmulatorDigis","Central"),
                                  l1UCTForwardJets        = cms.InputTag("uct2016EmulatorDigis","Forward"),
                                  stage2Taus              = cms.InputTag("l1extraParticles","Tau"),
                                  stage2IsoTaus           = cms.InputTag("l1extraParticles","IsoTau"),
                                  stage2DigisTaus         = cms.InputTag("caloStage2Digis", "Tau"),
                                  gmtStage2Digis          = cms.InputTag("simGmtStage2Digis"),
                                  genJets                 = cms.InputTag("slimmedGenJets"),
                                  isData                  = cms.bool(True),
                                  folderName              = cms.untracked.string("Stage3Regions")
)
#BXVector<l1t::Muon>                   "simGmtStage2Digis"         ""                "HLT"
