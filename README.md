# bbtoDijet
Modifying 2 CMS high-level triggers to lower the mass threshold while conserving the rate.

## Setup 
Use the [SWGuideGlobalHLT Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#Preparing_a_80X_CMSSW_developer/ "Preparing a working area for 80X") and [TSG tutorials](https://indico.cern.ch/event/520258/ "Trigger Tutorial") to prepare your working area area. The instructions below work on lxplus:

    cmsrel CMSSW_8_0_11
    cd CMSSW_8_0_11/src/
    cmsenv
    git cms-addpkg L1Trigger/L1TGlobal
    # HLT
    git cms-addpkg HLTrigger/Configuration
    git cms-merge-topic -u 14850
    git cms-merge-topic -u cms-tsg-storm:80XHLTJune15thTrain
    git cms-checkdeps -A -a
    ### Please make sure that the above compiles.
    git clone git@github.com:cms-steam/HLTrigger temp   
    cp -r temp/* HLTrigger/
    rm -fr temp/
    # L1 menu
    cp /afs/cern.ch/user/t/tmatsush/public/L1MenuDev/L1Menu_Collisions2016_v4.xml L1Trigger/L1TGlobal/data/Luminosity/startup/
    # compile -- may take a significant amount of time
    scram b -j 4
    cmsenv
    
Create a working directory:

    mkdir bbtoDijet
    cd bbtoDijet

Clone the repository:

    git clone https://github.com/avkhadiev/bbtoDijet.git

## Obtaining the trigger menu configuration file

Obtain the user menu from [ConfDB](https://cmsweb.cern.ch/confdb/ "HLT Configurations Explorer") at /users/aavkhadi/hlt_bTagDijet/V11. The procedure is similar to the one outlined in [HLTNtupleProductionSTEAM Twiki](https://twiki.cern.ch/twiki/bin/view/Sandbox/HLTNtupleProductionSTEAM#Create_CMSSW_config_files_user_m "Create CMSSW config files from a user menu"). 

    hltConfigFromDB --cff --configName /dev/CMSSW_8_0_0/GRun --nopaths --services -PrescaleService,-EvFDaqDirector,-FastMonitoringService > setup_cff.py
    hltGetConfiguration /users/aavkhadi/hlt_bbtoDijet/V11 --full --offline --data --unprescale --process TEST --l1Xml L1Menu_Collisions2016_v4.xml --globaltag 80X_dataRun2_HLT_v12 --input > hlt_bTagDijetV11.py
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

## What to do next

At this point, the following task are at hand:

1. Run the trigger on MC b-jet signal data to see what the distribution should look like. 
2. Estimate the fraction of b-jets in real data (requires summing data + signal, allow normalizations to float, and find best fit by maximizing the likelihood function)

## Current questions 

1. Do I get MC samples from David?
2. Once we know the shape of the true CSV distributution, what are the next steps?
3. How do we tell the trigger reached the target sensitivity? Do we look at the mass of two leading-pT jets?
4. To estimate the rate for our trigger, we will need to know the rate of our control trigger -- that requires running HLTPhysics jobs on control triggers to gauge their rate.
