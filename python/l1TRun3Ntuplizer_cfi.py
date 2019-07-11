import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("Run3Ntuplizer",
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  recoJets                = cms.InputTag("slimmedJets"),
                                  recoJetsAK8             = cms.InputTag("slimmedJetsAK8"),
                                  recoPtCut               = cms.double(10),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis'),
                                  l1UCTCentralJets        = cms.InputTag("uct2016EmulatorDigis","Central"),
                                  l1UCTForwardJets        = cms.InputTag("uct2016EmulatorDigis","Forward"),
                                  genJets                 = cms.InputTag("slimmedGenJets"),
                                  isData                  = cms.bool(True),
                                  folderName              = cms.untracked.string("Stage3Regions")
)
