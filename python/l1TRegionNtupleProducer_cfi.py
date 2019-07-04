import FWCore.ParameterSet.Config as cms


l1NtupleProducer = cms.EDAnalyzer("L1TRegionNtupleProducer",
                                  ecalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  hcalDigis               = cms.InputTag( 'l1tCaloLayer1Digis'),
                                  UCTRegion               = cms.InputTag('uct2016EmulatorDigis'),
                                  folderName              = cms.untracked.string("Stage3Regions"),
)
