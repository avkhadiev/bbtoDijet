# bbtoDijet
Modifying 2 CMS high-level triggers to lower the mass threshold while conserving the rate.

## Setup 
Use the [SWGuideGlobalHLT Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#Preparing_a_80X_CMSSW_developer/ "Preparing a working area for 80X") to prepare your working area area:

    cmsrel CMSSW_8_0_11
    cd CMSSW_8_0_11/src/
    cmsenv
    git cms-addpkg HLTrigger/Configuration
    git cms-merge-topic -u cms-tsg-storm:80XHLTJuly29thTrain
    git cms-checkdeps -A -a
    scram b -j 4
    cmsenv

Clone the repository:

    git clone https://github.com/avkhadiev/bbtoDijet.git

## Obtaining the trigger menu configuration file

Obtain the user menu from [ConfDB](https://cmsweb.cern.ch/confdb/ "HLT Configurations Explorer") at /users/aavkhadi/hlt_bbtoDijet/V11. The procedure is similar to the one outlined in [HLTNtupleProductionSTEAM Twiki](https://twiki.cern.ch/twiki/bin/view/Sandbox/HLTNtupleProductionSTEAM#Create_CMSSW_config_files_user_m "Create CMSSW config files from a user menu"). 

    hltConfigFromDB --cff --configName /dev/CMSSW_8_0_0/GRun --nopaths --services -PrescaleService,-EvFDaqDirector,-FastMonitoringService > setup_cff.py
    hltGetConfiguration /users/aavkhadi/hlt_bbtoDijet/V11 --full --offline --data --unprescale --process TEST --l1-emulator 'Full' --l1Xml L1Menu_Collisions2016_v4.xml --globaltag 80X_dataRun2_HLT_v12 --input > hlt_bbtoDijetV11.py
    # Edit the config file and add the following line just after 'process = cms.Process( "TEST" )': process.load("setup_cff")
