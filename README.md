# bbtoDijet
Modifying 2 CMS high-level triggers to lower the mass threshold while conserving the rate.

## Setup 
Use the [SWGuideGlobalHLT Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#Preparing_a_80X_CMSSW_developer/ "Preparing a working area for 80X") to prepare your working area area:

    cmsrel CMSSW_8_0_11
    cd CMSSW_8_0_11/src/
    cmsenv
    # the steps below take a significant amount of time
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
    hltGetConfiguration /users/aavkhadi/hlt_bbtoDijet/V11 --full --offline --data --unprescale --process TEST --l1Xml L1Menu_Collisions2016_v4.xml --globaltag 80X_dataRun2_HLT_v12 --input > hlt_bbtoDijetV11.py
    # Edit the config file and add the following line just after 'process = cms.Process( "TEST" )': process.load("setup_cff")
    
## Workflow

We modify the following triggers:

1. HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v2
2. HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v2

They can be found on [ConfDB](https://cmsweb.cern.ch/confdb/ "HLT Configurations Explorer") at /dev/CMSSW_8_0_0/GRun/V123. 
The first trigger has lower pT working point (at least two jets should have pT of 100 GeV), and higher CSV working point (at least two jets should have an online CSV discriminant of .86: 100 - 14 = 86, hence p014 in the name). The second trigger has a higher pT of 160 GeV and a lower CSV working point of .78 (typo in the name: should be p022, not p026).

For both triggers, we modify the following three parameters: pT, abs(eta), and CSV working point, all mainly in the relaxation direction. This is augmented by efficiency studies. We produce turn-on curves as functions of the three parameters, scanning around currently chosen values. We intend to lower the mass threshold (sum of masses of the two leading jets) below 750 GeV. 

For the efficiency studies, we use the following control triggers:

1. HLT_DoubleJetsC100_p014_DoublePFJetsC100MaxDeta1p6_v2 and HLT_DoubleJetsC100_p026_DoublePFJetsC160_v2 -- the triggers of interest with the CSV tagger removed, for b-tagging efficiency studies. 
2. HLT_Mu50_v3 as a control trigger for the kinematics efficiency studies.
