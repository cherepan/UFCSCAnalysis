from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'RAW2DIGI_RealVal'
config.General.workArea =  'crab_area_ZMMRelVal_20_08_2024_v2'
config.General.transferLogs = True


config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'makeLocalReco_TryRunIIConditions.py'
config.JobType.maxMemoryMB      = 8000
config.JobType.maxJobRuntimeMin = 2750
config.JobType.numCores = 4

config.section_("Data")

#config.Data.inputDataset = '/RelValZMM_14/cherepan-RelValZMM_CMSSW_13_3_0_FullTier_16_08_2024_mcRun3_2023-faac953493000e6dbd6a65acab9b9ebc/USER'
config.Data.inputDataset = '/RelValZMM_14/cherepan-RelValZMM_CMSSW_13_3_0_FullTier_19_08_2024_WinterValidation-faac953493000e6dbd6a65acab9b9ebc/USER'

config.Data.inputDBS  = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits  = 100
#config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/cherepan'
config.Data.outputDatasetTag = 'RelValZMM_CMSSW_13_3_0_FullTier_20_08_2024_AddUF_LR_v2'


config.section_("Site")
##config.Site.whitelist = ['T2_US_Wisconsin','T2_US_Purdue','T1_US_FNAL']
##config.Data.ignoreLocality = True
#config.Site.storageSite = 'T2_US_Florida'
config.Site.whitelist = ['T2_US_Florida']#,'T2_US_Wisconsin','T2_US_Purdue','T1_US_FNAL']
config.Data.ignoreLocality = True
config.Site.storageSite = 'T2_US_Florida'

