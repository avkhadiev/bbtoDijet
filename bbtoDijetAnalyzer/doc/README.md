Writing a New Analyzer

Example code:

 - [HLTAnalyzer.cc](https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/HLTrigger/HLTanalyzers/src/HLTAnalyzer.cc/ "HLT Analyzer")
 - [HLTBitAnalyzer.cc](https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/HLTrigger/HLTanalyzers/src/HLTBitAnalyzer.cc/ "HLT Bit Analyzer")
 - My SimpleHLTAnalyzer
 - [HLTScoutingCaloProducer.cc](https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/de/d43/HLTScoutingCaloProducer_8cc_source.html "HLT Scouting Calo Producer")

WTD:

 1. Add all necessary #include's
 2. Add a getCollection function similar to [HLTBitAnalyzer.cc](https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/HLTrigger/HLTanalyzers/src/HLTBitAnalyzer.cc/ "HLT Bit Analyzer")
 3. Write main method: define parameters, use TFileService, set up analysis
 4. Write the analyze method: define handles, getCollection, run the analysis, fill the tree
 5. Make sure all files are closed in the end to avoid excessive memory use
 6. Edit fillDescriptions like in [HLTScoutingCaloProducer.cc](https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/de/d43/HLTScoutingCaloProducer_8cc_source.html "HLT Scouting Calo Producer")
 
What to include in the ntuples:

1. Trigger results
2. CSV Online
3. CSV Offline 
4. Run, Lumi, Events
5. Online & Offline pt, eta, phi, mass
