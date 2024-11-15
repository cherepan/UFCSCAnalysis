######################################################################
#                                                                    #
# This version works with CMSSW_13_3_X   MC tested                   #
#                                                                    #
######################################################################
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()




doUnpacking = bool(True)

########## Options ############
runUF_localReco    = bool(True)    # True if run CSCSegAlgoUF; False - default

isDATA             = bool(False)
isRAW              = bool(True)
isDIGI             = bool(True)

isSIM              = bool(True)
isGEN              = bool(True)

isLocalRECO        = bool(True)
isFullRECO         = bool(True)

addMuonInfo        = bool(True)
addTrackInfo       = bool(False)
addRecHitInfo      = bool(True)
addSegmentInfo     = bool(True)
addTriggerInfo     = bool(True)

addDigiInfo        = bool(True)
addTimeMonitoringInfo = bool(True)
addCalibrationInfo    = bool(False)

maxEvents = -1


MCGlobalTag='133X_mcRun3_2024_realistic_v7'
DataGlobalTag='124X_dataRun3_PromptAnalysis_v1'




doDebug = bool(False)

###############################

### Debug Printing ###

if not isDATA:
    print('Sample Type: MC')
else:
    print("Sample Type: Data")
                



    
#####################
process = cms.Process("UFCSCRootMaker")

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
#process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.suppressWarning.append('classByHitsGlb') # kill stupid RPC hit associator warning

                                                       
######################### Frontier Conditions #########################
# Conditions data for calibration and alignment                       #
# are defined in the Offline Conditions Database (ORCOF),             #
# which is read in CMSSW applications via Frontier caching servers.   #
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions  #
#######################################################################


process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration/StandardSequences/MagneticField_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("UFCSCAnalysis.UFCSCRootMaker.cscRootMaker_cfi")
###
if not isDATA:
    process.GlobalTag.globaltag=MCGlobalTag 
    
else:
    process.GlobalTag.globaltag=DataGlobalTag


process.out = cms.OutputModule("PoolOutputModule",
                               #fileName = cms.untracked.string('test.root'),
                               # save only events passing the full path
                               outputCommands = cms.untracked.vstring('drop *')
                               )





process.cscRootMaker.cscRecHitTagSrc = cms.untracked.InputTag('csc2DRecHits','','localRecoUF')
process.cscRootMaker.cscSegTagSrc    = cms.untracked.InputTag('cscSegments','','localRecoUF')

if(not runUF_localReco):
    process.cscRootMaker.cscRecHitTagSrc = cms.untracked.InputTag('csc2DRecHits','','RECO')
    process.cscRootMaker.cscSegTagSrc    = cms.untracked.InputTag('cscSegments','','RECO')




outputFileName = 'CSC_UF_Ntuple_SegmentAlgoUF.root'
if(not runUF_localReco):
    outputFileName = 'CSC_UF_Ntuple_SegmentAlgoDefault.root'

    
process.TFileService = cms.Service("TFileService",

                                   fileName = cms.string("file:"+outputFileName)

                                   )

# Primary Vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                                  src = cms.InputTag('offlinePrimaryVertices'),
                                                  cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                                  filter = cms.bool(True)
                                                  )

##########  debug block #####
txt_file = open('Events_ToDebug.txt', 'r')
#if(not CleanChambers):
#    txt_file = open('Events_with_Noisy_CSC_with_muon.txt', 'r')

print("Getting events to process from:", txt_file)
event_list = [event.split('-')[0].strip() for event in txt_file.readlines()]
event_list= list(dict.fromkeys(event_list))
print(" --> running on %i events" % len(event_list) )
##########  debug block #####





process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             # processingMode = cms.untracked.string('RunsAndLumis'),
                             eventsToProcess     = cms.untracked.VEventRange(event_list),
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                             fileNames = cms.untracked.vstring(),      
                             )

if isDATA:
    process.source.fileNames = cms.untracked.vstring(
        '/store/user/cherepan/RelValSingleMuPt1/RelValZMM_14_CMSSW_12_4_0RAW2DIGI_L1Reco_RECO_RECOSIM_UF_SingleMu1_UF_LocalRECO/230908_124247/0000/SingleMu1Pt_14_UF_CSCSegmentBuilder_1.root',
        '/store/user/cherepan/RelValSingleMuPt1/RelValZMM_14_CMSSW_12_4_0RAW2DIGI_L1Reco_RECO_RECOSIM_UF_SingleMu1_UF_LocalRECO/230908_124247/0000/SingleMu1Pt_14_UF_CSCSegmentBuilder_2.root',
        '/store/user/cherepan/RelValSingleMuPt1/RelValZMM_14_CMSSW_12_4_0RAW2DIGI_L1Reco_RECO_RECOSIM_UF_SingleMu1_UF_LocalRECO/230908_124247/0000/SingleMu1Pt_14_UF_CSCSegmentBuilder_3.root'
)
else:
    process.source.fileNames = cms.untracked.vstring(


#  uncomment later
#        'file:/eos/user/c/cherepan/CSC/Test_AddUFLocalReco.root'
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_1.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_2.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_3.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_4.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_5.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_6.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_7.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_8.root',
        'file:/eos/user/c/cherepan/CSC/Test2/AddingUFLR/AddUFLocalReco_9.root'
    )
    
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True)
                                     #,SkipEvent = cms.untracked.vstring('ProductNotFound')
                                     )

### Physics Declared Filter (for data)
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

### No scraping
process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(doDebug),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.2)
                                  )


#process.LumiCorrectionSource = cms.ESSource("LumiCorrectionSource",
                                            #authpath=cms.untracked.string('/afs/cern.ch'file:/cms/lumi/DB'),
                                            #connect=cms.string('oracle:/'file:/cms_orcon_adg'file:/cms_lumi_prod')
#                                            connect=cms.string('frontier://LumiCalc_LUMI_PROD')
                                            #normtag=cms.untracked.string('HFV2a')
                                            #datatag=cms.untracked.string('v3')
#                                            )



#EventCount                                                                                                                             
process.nEventsTotal = cms.EDProducer("EventCountProducer")

#RM

process.cscRootMaker.isFullRECO = cms.untracked.bool(isFullRECO)
process.cscRootMaker.isLocalRECO = cms.untracked.bool(isLocalRECO)
process.cscRootMaker.isGEN = cms.untracked.bool(isGEN)
process.cscRootMaker.isSIM = cms.untracked.bool(isSIM)
process.cscRootMaker.isRAW = cms.untracked.bool(isRAW)
process.cscRootMaker.isDIGI = cms.untracked.bool(isDIGI)
process.cscRootMaker.isDATA = cms.untracked.bool(isDATA)
process.cscRootMaker.addMuons = cms.untracked.bool(addMuonInfo)
process.cscRootMaker.addTracks = cms.untracked.bool(addTrackInfo)
process.cscRootMaker.addRecHits = cms.untracked.bool(addRecHitInfo)
process.cscRootMaker.addSegments = cms.untracked.bool(addSegmentInfo)
process.cscRootMaker.addTrigger = cms.untracked.bool(addTriggerInfo)
process.cscRootMaker.addDigis = cms.untracked.bool(addDigiInfo)
process.cscRootMaker.addTimeMonitoring = cms.untracked.bool(addTimeMonitoringInfo)
process.cscRootMaker.addCalibrations = cms.untracked.bool(addCalibrationInfo)

process.unpack = cms.Sequence(process.muonCSCDigis*process.gtDigis)

process.p = cms.Path(
    process.nEventsTotal*process.cscRootMaker
)

if not isDATA:
    process.p.replace(process.nEventsTotal,process.nEventsTotal)

if isDATA:
#    process.load("CondCore.CondDB.CondDB_cfi")

    process.p.replace(process.nEventsTotal,process.nEventsTotal)

if doUnpacking:
    process.p.replace(process.nEventsTotal,process.unpack*process.nEventsTotal)

