############################
#                          #
#    JetHT Run 2016B       #
#                          #
############################

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

name = 'eff_bTagDijetV11' # will be part of the work area name and the storage subdir name
runNom = '274968-274969,274998,274999' #

# For information on config parameters, see
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

# section General
config.General.requestName = name
config.General.workArea = 'crab_test_' + name
config.General.transferOutputs = True
config.General.transferLogs = True

# section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hlt_bTagDijetV11.py'
config.JobType.outputFiles = ['hlt_bTagDijetV11.root']
config.JobType.numCores = 16
# section Data
config.Data.inputDataset = '/JetHT/Run2016B-PromptReco-v2/AOD' # '/HLTPhysics/Run2016B-PromptReco-v2/AOD'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 3  # use crab submit --dryrun *.py to find optimal splitting
config.Data.secondaryInputDataset = '/JetHT/Run2016B-v2/RAW'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt' # specifes good lumi sections to be used
config.Data.totalUnits = -1 # analyze all events after applying the lumi mask
config.Data.runRange = runNom
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False # no need to publish the results
config.Data.outputDatasetTag = name
config.Data.ignoreLocality = True

# section Site
config.Site.storageSite = 'T3_US_FNALLPC'
