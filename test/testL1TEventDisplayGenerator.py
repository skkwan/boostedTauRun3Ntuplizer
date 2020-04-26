import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest1")

#import EventFilter.L1TXRawToDigi.util as util

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.parseArguments()

#def formatLumis(lumistring, run) :
#    lumis = (lrange.split('-') for lrange in lumistring.split(','))
#    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
#    return ['-'.join(l) for l in runlumis]


#if len(options.inputFiles) is 0 and options.inputFileList is '' :
#    inputFiles = util.getFilesForRun(options.runNumber, options.dataStream)
#elif len(options.inputFileList) > 0 :
#    with open(options.inputFileList) as f :
#        inputFiles = list((line.strip() for line in f))
#else :
#    inputFiles = cms.untracked.vstring(options.inputFiles)
#if len(inputFiles) is 0 :
#    raise Exception('No files found for dataset %s run %d' % (options.dataStream, options.runNumber))
#print 'Ok, time to analyze'
#
#secondaryMap = {
#    "root://xrootd-cms.infn.it//store/data/Run2015D/SingleMuon/MINIAOD/05Oct2015-v1/10000/04EDCDA3-916F-E511-AD11-0025905938B4.root": [
#        'root://xrootd-cms.infn.it//store/data/Run2015D/SingleMuon/RAW/v1/000/257/487/00000/28774BB9-EF66-E511-A328-02163E011CC3.root',
#        'root://xrootd-cms.infn.it//store/data/Run2015D/SingleMuon/RAW/v1/000/257/487/00000/8A42A2DA-EF66-E511-AC62-02163E012988.root',
#        'root://xrootd-cms.infn.it//store/data/Run2015D/SingleMuon/RAW/v1/000/257/487/00000/C24B52CB-EF66-E511-969B-02163E0119F6.root']
#}

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

## To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_2016_v2_1_cfi')
## To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

#process.load("L1Trigger.Run3Ntuplizer.l1BoostedJetStudies_cfi")

##if options.data is False :
#process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
#process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
#        cms.InputTag('hcalDigis'),
#        cms.InputTag('hcalDigis')
#)
#process.L1TReEmul = cms.Sequence(process.simHcalTriggerPrimitiveDigis * process.SimL1Emulator)

process.load("L1Trigger.Run3Ntuplizer.l1TEventDisplayGenerator_cfi")
process.l1NtupleProducer.ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives","L1TCaloSummaryTest1")
process.l1NtupleProducer.hcalDigis = cms.InputTag("hcalDigis")

process.uct2016EmulatorDigis.useECALLUT = cms.bool(True)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(True)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/SUSYGluGluToBBHToBB_NarrowWidth_M-1200_TuneCP5_13TeV-pythia8/MINIAODSIM/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/110000/1B0FC5E0-ACF7-0C49-8262-1F95F23C896B.root'
),
                            secondaryFileNames = cms.untracked.vstring(
                                'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18DR/SUSYGluGluToBBHToBB_NarrowWidth_M-1200_TuneCP5_13TeV-pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/110000/4187DAD9-3092-E448-9F4E-6B3616479547.root',
                                'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18DR/SUSYGluGluToBBHToBB_NarrowWidth_M-1200_TuneCP5_13TeV-pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/110000/2FE6F8A3-9F3C-E14B-8F5F-F45BBE34AA30.root',
                                'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18DR/SUSYGluGluToBBHToBB_NarrowWidth_M-1200_TuneCP5_13TeV-pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/110000/44D65997-B07E-8143-8B7D-B15A04DDC388.root',
                                'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18DR/SUSYGluGluToBBHToBB_NarrowWidth_M-1200_TuneCP5_13TeV-pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/110000/6B447B66-341D-E742-A73A-AB39EEDD25CA.root',
                                'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18DR/SUSYGluGluToBBHToBB_NarrowWidth_M-1200_TuneCP5_13TeV-pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/110000/8A7F5240-F850-934B-A3AE-A0FFC746FDCB.root',
                            )
)
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:734","1:961","1:966","1:982","1:966")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TDisplayEvent.root")
)

#if options.data is True :
#    process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)
#else :

process.p = cms.Path(process.RawToDigi*process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("dumpout.root"),
#    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
#)
#process.e = cms.EndPath(process.out)
