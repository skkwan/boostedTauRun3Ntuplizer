import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")

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

def formatLumis(lumistring, run) :
    lumis = (lrange.split('-') for lrange in lumistring.split(','))
    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
    return ['-'.join(l) for l in runlumis]

print 'Getting files for run %d...' % options.runNumber
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
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1TRun3Ntuplizer_cfi")

process.l1NtupleProducer.isData = cms.bool(False)

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring(inputFiles)#,
                            #secondaryFileNames = cms.untracked.vstring(secondaryMap[options.inputFiles[0]])
                            fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/ABCCF610-63D3-7046-AA37-4D8B3F1D3494.root'
),
                            secondaryFileNames = cms.untracked.vstring(
                                '/store/mc/RunIIAutumn18DR/VBFHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/018EB922-D39E-C045-B7B9-F07C91216849.root',
                                '/store/mc/RunIIAutumn18DR/VBFHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/018EB922-D39E-C045-B7B9-F07C91216849.root',
                                '/store/mc/RunIIAutumn18DR/VBFHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/3FFEF5C8-A16B-CB46-98A9-49B8674E96E7.root',
                                '/store/mc/RunIIAutumn18DR/VBFHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/90091689-95C4-E44A-B681-2AFCE2F2F0DB.root',
                                '/store/mc/RunIIAutumn18DR/VBFHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/BC74B35D-B1CB-214B-9114-ECFA365539B9.root',
                                '/store/mc/RunIIAutumn18DR/VBFHToTauTau_M125_13TeV_powheg_pythia8/GEN-SIM-RAW/FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/60000/F25C6600-C8C9-BB46-8E33-7F755B79A93B.root'
                            )
)

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:3375","1:1443","1:2960","1:3805","1:3381","1:229","1:3531","1:583","1:1120","1:504","1:1454","1:874","1:1078","1:1087","1:121","1:1448","1:1456","1:177","1:2964")


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent.root"),
    outputCommands = cms.untracked.vstring('keep *')
    #outputCommands = cms.untracked.vstring('drop *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-VBF.root")
)

process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

process.e = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.p,process.e)
process.schedule = cms.Schedule(process.p)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

dump_file = open('dump.py','w')
dump_file.write(process.dumpPython())
