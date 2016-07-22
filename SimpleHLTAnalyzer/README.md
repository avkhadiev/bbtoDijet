# SimpleHLTAnalyzer 

Simple HLT Analyzer is an ED analyzer the reference to which is included in the configuration file for a trigger menu. 
The analyzer creates an ouput tree, which branches contain the necessary information about each event a job was run over. 
The branches include on- and offline CSV discriminants per jet, run numbers, luminosity blocks, and jet pTs, etas, phis, and masses. 

## Setup

You need to 

    scram b -j 4 

to compile the analyzer code locally. You need to compile this way each time you make changes to the analyzer code.
The code itself is stored in macros/. In test/, you can find trigger menu configuration files, and 
corresponding configuration files for CRAB jobs. In macros/ you fill find a Python script that makes efficiency turn-on 
curves, provided a .root file with branches -- the output of the ED analyzer. 

## Warning

The analyzer does not use memory sparingly, resulting in ~ 10% jobs fail due to excessive memory use. An attempt to write a new analyzer from scratch is being made, see bbtoDijetAnalyzer.
