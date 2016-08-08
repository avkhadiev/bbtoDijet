############################
#                          #
#   HLTPhysics Run 2016B   #
#                          #
############################

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

name = 'bbtoDijetV11_HLT' # will be part of the work area name and the storage subdir name
runNom = '274998' #

# For information on config parameters, see
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

# section General
config.General.requestName = name
config.General.workArea = 'crab_test_' + name
config.General.transferOutputs = True
config.General.transferLogs = True

# section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hlt_bTagDijetV11_HLT.py'
config.JobType.outputFiles = ['hlt_bTagDijetV11_HLT.root']
config.JobType.numCores = 8
# section Data
config.Data.inputDataset = '/HLTPhysics/Run2016B-PromptReco-v2/AOD'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 3  # use crab submit --dryrun *.py to find optimal splitting
config.Data.secondaryInputDataset = '/HLTPhysics/Run2016B-v2/RAW'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt' # specifes good lumi sections to be used
config.Data.totalUnits = -1 # analyze all events after applying the lumi mask
config.Data.runRange = runNom
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False # no need to publish the results
config.Data.outputDatasetTag = name
config.Data.ignoreLocality = True

# section Site
config.Site.storageSite = 'T3_US_FNALLPC'
