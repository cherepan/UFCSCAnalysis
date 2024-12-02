## Dump  10  events in CSC rechit builder - Tim Cox - 07.11.2012
## This version runs in 6_0_1_PostLS1 on a simulated data DIGI sample.

import FWCore.ParameterSet.Config as cms

process = cms.Process("localRecoUF")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load('Configuration.StandardSequences.Services_cff')

##################################################
# --- MATCH GT TO RELEASE AND DATA SAMPLE
#process.GlobalTag.globaltag = "POSTLS161_V11::All"
#process.GlobalTag.globaltag = "106X_dataRun2_v32"


#process.GlobalTag.globaltag = "124X_dataRun3_PromptAnalysis_v1"
#process.GlobalTag.globaltag = "124X_mcRun3_2022_realistic_v12"
process.GlobalTag.globaltag = "133X_mcRun3_2024_realistic_v7"


#process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v15_L1v1"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet(  500 )
isSIMDIGI = bool(False)
isRAW = bool(True)
###################################################


CleanChambers = True;




txt_file = open('Events_with_Clean_CSC_with_muon.txt', 'r')


#txt_file = open('EventsAndChambersWithTwoSegments.txt', 'r')

#txt_file = open('EventsAndChambersWithFourSegments.txt', 'r')
#txt_file = open('Events_with_more_that_one_segment.txt', 'r')
#txt_file = open('Events_ToDebug.txt', 'r')


if(not CleanChambers):
    txt_file = open('Events_with_Noisy_CSC_with_muon.txt', 'r')

print("Getting events to process from:", txt_file)
event_list = [event.split('-')[0].strip() for event in txt_file.readlines()]
event_list= list(dict.fromkeys(event_list))
print(" --> running on %i events" % len(event_list), ' CleanChambers = ', CleanChambers)





process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring("ProductNotFound") )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source    = cms.Source("PoolSource",
                               noEventSort         = cms.untracked.bool(True),
                               eventsToProcess     = cms.untracked.VEventRange(event_list),
                               duplicateCheckMode  = cms.untracked.string('noDuplicateCheck'),
                               fileNames = cms.untracked.vstring(
#                                    'file:/eos/user/c/cherepan/CSC/InputFilesZMMRelVal_13_3_0/a95425f8-0854-4252-bc17-42ccb35ab3fd.root'

                                   'file:/eos/user/c/cherepan/CSC/Test2/ZMMRelVal_13_3_0_RecoTier_1.root'
#                                   'file:/eos/user/c/cherepan/CSC/InputFilesZMMRelVal_13_3_0/744eef22-8474-4a45-9e08-01835b90cf62.root'
#                                   'file:/eos/user/c/cherepan/CSC/InputFilesZMMRelVal_13_3_0/1f317388-ad1b-463e-b053-5a69f4ad8364.root'

#                                    'file:/eos/user/c/cherepan/CSC/ZMMRelVal_13_3_0_FullTier_a95425f8-0854-4252-bc17-42ccb35ab3fd.root'


    )
)

# ME1/1A is  u n g a n g e d  Post-LS1
process.CSCGeometryESModule.useGangedStripsInME1a = True
##process.CSCGeometryESModule.debugV = True
##process.idealForDigiCSCGeometry.useGangedStripsInME1a = False

# Turn off some flags for CSCRecHitD that are turned ON in default config
process.csc2DRecHits.readBadChannels = cms.bool(False)
process.csc2DRecHits.readBadChannels = cms.bool(False)
process.csc2DRecHits.CSCUseTimingCorrections = cms.bool(False)
process.cscSegments.CSCUseTimingCorrections = cms.bool(False)
#process.csc2DRecHits.CSCUseTimingCorrections = cms.bool(False)
#process.csc2DRecHits.CSCUseGasGainCorrection = cms.bool(False)



# Switch input for CSCRecHitD to  s i m u l a t e d  digis  # i assume it is for MC only

process.csc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
process.csc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")

process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(1)


process.out = cms.OutputModule("PoolOutputModule",
                               fastCloning = cms.untracked.bool(False),
#                               fileName = cms.untracked.string('/eos/user/c/cherepan/CSC/UF_tuples/SingleMuon_RAW-RECO_ZMu-12Nov2019_UL2018_CSCSegmentBuilder_UF_testRun.root'),
#                               fileName = cms.untracked.string('RelValZMM_14_RU_CSCSegmentBuilder.root'),
#                               fileName = cms.untracked.string('/eos/user/c/cherepan/CSC/LocalReco_output_UF/Single_14_UF_CSCSegmentBuilder.root'),
#                               fileName = cms.untracked.string('/eos/user/c/cherepan/CSC/Test_AddUFLocalReco.root'),
#                               fileName = cms.untracked.string('/eos/user/c/cherepan/CSC/AddUFLocalReco_a95425f8-0854-4252-bc17-42ccb35ab3fd.root'),
                               fileName = cms.untracked.string('/eos/user/c/cherepan/CSC/AddUFLocalReco_1.root'),


#                               fileName = cms.untracked.string('DY_MUMURelVal_UF_CSCSegmentBuilder.root'),
#                               fileName = cms.untracked.string('DYJetsToLL_M-50_RU_CSCSegmentBuilder.root'),
                               outputCommands = cms.untracked.vstring('keep *')
                               )


# --- TO ACTIVATE LogTrace IN CSCRecHitD NEED TO COMPILE IT WITH scram b -j8 USER_CXXFLAGS="-DEDM_ML_DEBUG"
# LogTrace output goes to cout; all other output to "junk.log"
 
from Configuration.DataProcessing.RecoTLR import customiseDataRun2Common

#call to customisation function customiseDataRun2Common imported from Configuration.DataProcessing.RecoTLR
process = customiseDataRun2Common(process)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
# module label is something like "muonCSCDigis"...

#process.MessageLogger.debugModules = cms.untracked.vstring("*")

#process.MessageLogger.debugModules = ["CSCRecHit", "CSCRecHitDBuilder","CSCMake2DRecHit","CSCHitFromStripOnly"]
process.MessageLogger.debugModules = [ "CSCSegAlgoUF"]

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5)
#process.MessageLogger.destinations = cms.untracked.vstring("cout","junk")
process.MessageLogger.cerr.threshold = "DEBUG"
process.MessageLogger.cerr.enable = True

process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string("DEBUG"),
    default   = cms.untracked.PSet( limit = cms.untracked.int32(1)  ),
    FwkReport = cms.untracked.PSet( limit = cms.untracked.int32(5) )
##    , CSCGeometry = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#    , CSCRecHit = cms.untracked.PSet( limit = cms.untracked.int32(20) )
    , CSCRecHitDBuilder = cms.untracked.PSet( limit = cms.untracked.int32(20) )
#    , CSCMake2DRecHit = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#    , CSCHitFromStripOnly = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
##    , CSCRecoConditions = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)




process.MessageLogger = cms.Service("MessageLogger",
     destinations  = cms.untracked.vstring(
         'detailedInfo'
         ,'critical'
         ,'cerr'
     )
)


process.load("SimMuon.CSCDigitizer.muonCSCDigis_cfi")
# Path and EndPath def
process.unpack = cms.Path(process.muonCSCDigis * process.gtDigis)
process.reco = cms.Path(process.csc2DRecHits * process.cscSegments )




#process.reco = cms.Path(process.cscSegments )
#process.reco = cms.Path(process.reconstruction)
process.out_step = cms.EndPath(process.out)

# Schedule definition
process.schedule = cms.Schedule(process.reco, process.out_step)


if isSIMDIGI:
    process.reco.replace(process.csc2DRecHits, process.simMuonCSCDigis* process.gtDigis * process.csc2DRecHits * process.cscSegments )

if isRAW:
#    process.reco.replace(process.csc2DRecHits, process.muonCSCDigis * process.gtDigis * process.csc2DRecHits * process.cscSegments)
    process.reco.replace(process.csc2DRecHits, process.muonCSCDigis * process.gtDigis * process.csc2DRecHits * process.cscSegments)
