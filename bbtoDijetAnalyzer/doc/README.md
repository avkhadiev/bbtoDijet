Writing a New Analyzer

Example code:

 - [HLTAnalyzer.cc](https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/HLTrigger/HLTanalyzers/src/HLTAnalyzer.cc/ "HLT Analyzer")
 - [HLTBitAnalyzer.cc](https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/HLTrigger/HLTanalyzers/src/HLTBitAnalyzer.cc/ "HLT Bit Analyzer")
 - My SimpleHLTAnalyzer
 - [HLTScoutingCaloProducer.cc](https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/de/d43/HLTScoutingCaloProducer_8cc_source.html "HLT Scouting Calo Producer")
 - Also, see other code in HLTrigger/HLTanalyzers directory on doxygen -- in particular, HLTJets.cc and HLTJets.h, and HLTBJet.cc and HLTBJet.h.

WTD:

- Running with CRAB on a small-ish sample to test memory management showed the excessive memory usage is not resolved. 
- Try reducing number of lumi sections per job
