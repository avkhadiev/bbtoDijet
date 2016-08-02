// -*- C++ -*-
//
// Package:    bbtoDijet/bbtoDijetAnalyzer
// Class:      bbtoDijetAnalyzer
// 
/**\class bbtoDijetAnalyzer bbtoDijetAnalyzer.cc bbtoDijet/bbtoDijetAnalyzer/plugins/bbtoDijetAnalyzer.cc

 Description: collect b-tagging data (online and offline), calo and pf jet kinematics, and trigger results 

 Implementation:
    - set up tree variables during construction
    - clear (memset) variables before each event
    - collect data in analyze() per event
    - delete tree variables during destruction
    FIXME: branches for trigger results are created in analyze();
        corresponding arrays are not cleared nor deleted. 
        This can lead to excessive memory use on some CRAB jobs, need testing to verify. 
    FIXME: getting data for calo and pf jets, and getting data for online and offline CSV tags is very similar; 
        writing two functions to get the event information will simplify code
*/
//
// Original Author:  Artur Avkhadiev
//         Created:  Fri, 22 Jul 2016 08:05:25 GMT
//
//

// system include file
#include <memory>                                                     // included by default
#include <TTree.h>                                                    // added in
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"                  // included by default
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"                    // added in
#include "FWCore/Framework/interface/ESHandle.h"                      // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/dd/d5c/HLTBitAnalyzer_8h_source.html
#include "FWCore/Framework/interface/ConsumesCollector.h"             // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/d04/HLTJets_8h_source.html
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/df/d03/HLTEventAnalyzerRAW_8cc_source.html
#include "FWCore/MessageLogger/interface/MessageLogger.h"             // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d5/d74/HLTBitAnalyzer_8cc_source.html
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTanalyzers/interface/EventHeader.h"             // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/dd/d5c/HLTBitAnalyzer_8h_source.html
#include "HLTrigger/HLTanalyzers/src/EventHeader.cc"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"            // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/dc0/HLTEventAnalyzerAOD_8cc_source.html
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/Handle.h"
// event content access
#include "DataFormats/Common/interface/TriggerResults.h"              // trigger-related information
#include "DataFormats/HLTReco/interface/TriggerEvent.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h"          // jet-related information 
#include "DataFormats/JetReco/interface/PFJetCollection.h"            // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/d04/HLTJets_8h_source.html
#include "DataFormats/JetReco/interface/Jet.h"                        // b-jet-tagging information
#include "DataFormats/BTauReco/interface/JetTag.h"                    // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/d2e/HLTBJet_8h_source.html  

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class bbtoDijetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit bbtoDijetAnalyzer(const edm::ParameterSet&);
      ~bbtoDijetAnalyzer();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void clear();
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
 
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) ;
      virtual void endRun  (edm::Run const&, edm::EventSetup const&) ;

     // ----------member data ---------------------------
      
     HLTConfigProvider                                       hltConfig_                    ;
     EventHeader                                             eventHeader_                  ; 
     std::string                                             processName_                  ;
     
     // input tags
     const edm::InputTag                                     triggerResultsTag_            ;
     const edm::InputTag                                     triggerEventTag_              ;
     const edm::InputTag                                     caloJetsTag_                  ;
     const edm::InputTag                                     pfJetsTag_                    ; 
     const edm::InputTag                                     bTagCSVOnlineTag_             ; 
     const edm::InputTag                                     bTagCSVOfflineTag_            ;
     
     // tokens                                              
     const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken_          ;
     const edm::EDGetTokenT<trigger::TriggerEvent>           triggerEventToken_            ;
     const edm::EDGetTokenT<reco::CaloJetCollection>         caloJetsToken_                ;
     const edm::EDGetTokenT<reco::PFJetCollection>           pfJetsToken_                  ;
     const edm::EDGetTokenT<reco::JetTagCollection>          bTagCSVOnlineToken_           ; 
     const edm::EDGetTokenT<reco::JetTagCollection>          bTagCSVOfflineToken_          ;
     
     // tree variables
     // array sizes
     static const size_t kMaxTriggerPass_     = 10000   ;
     static const size_t kMaxCaloJet_         = 10000   ;
     static const size_t kMaxPFJet_           = 10000   ;
     static const size_t kMaxBTagCSVOnline_   = 1000    ;
     static const size_t kMaxBTagCSVOffline_  = 1000    ;
     
     // set of variables for 
     // event-header information: included in eventHeader_
     // trigger-related information
     bool          firstEvent_                 ;   // will create new branch for trigger path if firstEvent_ = true
     unsigned int *triggerPass_                ;
     // jet-related information
     // CaloJet
     int           caloJetCounter_             ;     
     double       *caloJetPt_                  ;  
     float        *caloJetEta_                 ; 
     float        *caloJetPhi_                 ;  
     double       *caloJetEt_                  ;  
     float        *caloJetMass_                ;
     double       *caloJetEnergy_              ;
     double       *caloJetPx_                  ;
     double       *caloJetPy_                  ;
     double       *caloJetPz_                  ;
     // PFJet
     int           pfJetCounter_               ;     
     double       *pfJetPt_                    ;  
     float        *pfJetEta_                   ; 
     float        *pfJetPhi_                   ;  
     double       *pfJetEt_                    ;  
     double       *pfJetMass_                  ;
     double       *pfJetEnergy_                ;
     double       *pfJetPx_                    ;
     double       *pfJetPy_                    ;
     double       *pfJetPz_                    ;
     // b-jet-tagging information
     // bTagCSVOnline
     int           bTagCSVOnlineJetCounter_    ;     
     double       *bTagCSVOnlineJetPt_         ;
     float        *bTagCSVOnlineJetEta_        ;
     float        *bTagCSVOnlineJetPhi_        ;
     double       *bTagCSVOnlineJetEt_         ;
     double       *bTagCSVOnlineJetMass_       ;
     float        *bTagCSVOnline_              ;
     // bTagCSVOffline
     int           bTagCSVOfflineJetCounter_   ;
     double       *bTagCSVOfflineJetPt_        ;
     float        *bTagCSVOfflineJetEta_       ;
     float        *bTagCSVOfflineJetPhi_       ;
     double       *bTagCSVOfflineJetEt_        ;
     double       *bTagCSVOfflineJetMass_      ;
     float        *bTagCSVOffline_             ;
     
     TTree* outTree_;

};

//
// constructors and destructor
//
bbtoDijetAnalyzer::bbtoDijetAnalyzer(const edm::ParameterSet& iConfig) : 
  processName_                       (iConfig.getParameter<std::string>("processName")),  
  // initialize input tags 
  // (strings should match variable names in the *.py config file)
  triggerResultsTag_                 (iConfig.getParameter<edm::InputTag>("triggerResults")), 
  triggerEventTag_                   (iConfig.getParameter<edm::InputTag>("triggerEvent")),                     
  caloJetsTag_                       (iConfig.getParameter<edm::InputTag>("caloJets")),            
  pfJetsTag_                         (iConfig.getParameter<edm::InputTag>("pfJets")),             
  bTagCSVOnlineTag_                  (iConfig.getParameter<edm::InputTag>("bTagCSVOnline")),            
  bTagCSVOfflineTag_                 (iConfig.getParameter<edm::InputTag>("bTagCSVOffline")),            
  // initialize tokens
  triggerResultsToken_               (consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerEventToken_                 (consumes<trigger::TriggerEvent>(triggerEventTag_)), 
  caloJetsToken_                     (consumes<reco::CaloJetCollection>(caloJetsTag_)),
  pfJetsToken_                       (consumes<reco::PFJetCollection>(pfJetsTag_)),
  bTagCSVOnlineToken_                (consumes<reco::JetTagCollection>(bTagCSVOnlineTag_)),
  bTagCSVOfflineToken_               (consumes<reco::JetTagCollection>(bTagCSVOfflineTag_)), 
  // initialize tree variables
  firstEvent_                        (true)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
  
   /*  Setup the analysis to put the branch-variables into the tree. */
    
   // set of variables for
   // event-header information: included in eventHeader_ 
   // trigger-related information
   triggerPass_               = new unsigned int[kMaxTriggerPass_]  ;
   // jet-related information
   // CaloJet
   caloJetCounter_            = 0                                   ;
   caloJetPt_                 = new double[kMaxCaloJet_]            ;        
   caloJetEta_                = new float[kMaxCaloJet_]             ; 
   caloJetPhi_                = new float[kMaxCaloJet_]             ; 
   caloJetEt_                 = new double[kMaxCaloJet_]            ; 
   caloJetMass_               = new float[kMaxCaloJet_]             ;
   caloJetEnergy_             = new double[kMaxCaloJet_]            ;
   caloJetPx_                 = new double[kMaxCaloJet_]            ; 
   caloJetPy_                 = new double[kMaxCaloJet_]            ;
   caloJetPz_                 = new double[kMaxCaloJet_]            ;
   // PFJet
   pfJetCounter_              = 0                                   ;
   pfJetPt_                   = new double[kMaxPFJet_]              ;         
   pfJetEta_                  = new float[kMaxPFJet_]               ;  
   pfJetPhi_                  = new float[kMaxPFJet_]               ;  
   pfJetEt_                   = new double[kMaxPFJet_]              ;  
   pfJetMass_                 = new double[kMaxPFJet_]              ; 
   pfJetEnergy_               = new double[kMaxPFJet_]              ;
   pfJetPx_                   = new double[kMaxPFJet_]              ; 
   pfJetPy_                   = new double[kMaxPFJet_]              ;
   pfJetPz_                   = new double[kMaxPFJet_]              ;
   // b-jet-tagging information
   // bTagCSVOnline
   bTagCSVOnlineJetCounter_   = 0                                   ;
   bTagCSVOnlineJetPt_        = new double[kMaxBTagCSVOnline_]      ;
   bTagCSVOnlineJetEta_       = new float[kMaxBTagCSVOnline_]       ;
   bTagCSVOnlineJetPhi_       = new float[kMaxBTagCSVOnline_]       ;
   bTagCSVOnlineJetEt_        = new double[kMaxBTagCSVOnline_]      ;
   bTagCSVOnlineJetMass_      = new double[kMaxBTagCSVOnline_]      ;
   bTagCSVOnline_             = new float[kMaxBTagCSVOnline_]       ;
   // bTagCSVOffline 
   bTagCSVOfflineJetCounter_  = 0                                   ;  
   bTagCSVOfflineJetPt_       = new double[kMaxBTagCSVOffline_]     ;
   bTagCSVOfflineJetEta_      = new float[kMaxBTagCSVOffline_]      ;
   bTagCSVOfflineJetPhi_      = new float[kMaxBTagCSVOffline_]      ;
   bTagCSVOfflineJetEt_       = new double[kMaxBTagCSVOffline_]     ;
   bTagCSVOfflineJetMass_     = new double[kMaxBTagCSVOffline_]     ;
   bTagCSVOffline_            = new float[kMaxBTagCSVOffline_]      ;

   // open the tree file and initialize the tree
   edm::Service<TFileService> fs                                    ;
   outTree_ = fs->make<TTree>("efficiencyTree", "")                 ;

   // setup event header and HLT results analysis
   eventHeader_.setup(consumesCollector(), outTree_);

   // create branches
   // slash letter is ROOT's way for figuring out types:
   // I is int
   // F is float
   // D is double
   // sets of variables for
   // event-header information: included in eventHeader_ 
   // trigger-related information: created in ``analyze'' 
   // jet-related information
   // CaloJet
   outTree_->Branch("caloJetCounter",           &caloJetCounter_,            "caloJetCounter/I")                                   ;
   outTree_->Branch("caloJetPt",                 caloJetPt_,                 "caloJetPT[caloJetCounter]/D")                        ;    
   outTree_->Branch("caloJetEta",                caloJetEta_,                "caloJetEta[caloJetCounter]/F")                       ;    
   outTree_->Branch("caloJetPhi",                caloJetPhi_,                "caloJetPhi[caloJetCounter]/F")                       ;    
   outTree_->Branch("caloJetEt",                 caloJetEt_,                 "caloJetEt[caloJetCounter]/D")                        ;    
   outTree_->Branch("caloJetMass",               caloJetMass_,               "caloJetMass[caloJetCounter]/F")                      ;    
   outTree_->Branch("caloJetEnergy",             caloJetEnergy_,             "caloJetEnergy[caloJetCounter]/D")                    ;
   outTree_->Branch("caloJetPx",                 caloJetPx_,                 "caloJetPx[caloJetCounter]/D")                        ;
   outTree_->Branch("caloJetPy",                 caloJetPy_,                 "caloJetPy[caloJetCounter]/D")                        ;
   outTree_->Branch("caloJetPz",                 caloJetPz_,                 "caloJetPz[caloJetCounter]/D")                        ;
   // PFJet
   outTree_->Branch("pfJetCounter",             &pfJetCounter_,              "pfJetCounter/I")                                     ;
   outTree_->Branch("pfJetPt",                   pfJetPt_,                   "pfJetPT[pfJetCounter]/D")                            ;    
   outTree_->Branch("pfJetEta",                  pfJetEta_,                  "pfJetEta[pfJetCounter]/F")                           ;    
   outTree_->Branch("pfJetPhi",                  pfJetPhi_,                  "pfJetPhi[pfJetCounter]/F")                           ;    
   outTree_->Branch("pfJetEt",                   pfJetEt_,                   "pfJetEt[pfJetCounter]/D")                            ;    
   outTree_->Branch("pfJetMass",                 pfJetMass_,                 "pfJetMass[pfJetCounter]/D")                          ;    
   outTree_->Branch("pfJetEnergy",               pfJetEnergy_,               "pfJetEnergy[pfJetCounter]/D")                        ;
   outTree_->Branch("pfJetPx",                   pfJetPx_,                   "pfJetPx[pfJetCounter]/D")                            ;
   outTree_->Branch("pfJetPy",                   pfJetPy_,                   "pfJetPy[pfJetCounter]/D")                            ;
   outTree_->Branch("pfJetPz",                   pfJetPz_,                   "pfJetPz[pfJetCounter]/D")                            ;
   // b-jet-tagging information
   // bTagCSVOnline
   outTree_->Branch("bTagCSVOnlineJetCounter",  &bTagCSVOnlineJetCounter_,   "bTagCSVOnlineJetCounter/I")                          ;
   outTree_->Branch("bTagCSVOnlineJetPt",        bTagCSVOnlineJetPt_,        "bTagCSVOnlineJetPt[bTagCSVOnlineJetCounter]/D")      ;
   outTree_->Branch("bTagCSVOnlineJetEta",       bTagCSVOnlineJetEta_,       "bTagCSVOnlineJetEta[bTagCSVOnlineJetCounter]/F")     ;
   outTree_->Branch("bTagCSVOnlineJetPhi",       bTagCSVOnlineJetPhi_,       "bTagCSVOnlineJetPhi[bTagCSVOnlineJetCounter]/F")     ;
   outTree_->Branch("bTagCSVOnlineJetEt",        bTagCSVOnlineJetEt_,        "bTagCSVOnlineJetEt[bTagCSVOnlineJetCounter]/D")      ;
   outTree_->Branch("bTagCSVOnlineJetMass",      bTagCSVOnlineJetMass_,      "bTagCSVOnlineJetMass[bTagCSVOnlineJetCounter]/D")    ;
   outTree_->Branch("bTagCSVOnline",             bTagCSVOnline_,             "bTagCSVOnline[bTagCSVOnlineJetCounter]/F")           ;
   // bTagCSVOffline
   outTree_->Branch("bTagCSVOfflineJetCounter", &bTagCSVOfflineJetCounter_,  "bTagCSVOfflineJetCounter/I")                         ;
   outTree_->Branch("bTagCSVOfflineJetPt",       bTagCSVOfflineJetPt_,       "bTagCSVOfflineJetPt[bTagCSVOfflineJetCounter]/D")    ;
   outTree_->Branch("bTagCSVOfflineJetEta",      bTagCSVOfflineJetEta_,      "bTagCSVOfflineJetEta[bTagCSVOfflineJetCounter]/F")   ;
   outTree_->Branch("bTagCSVOfflineJetPhi",      bTagCSVOfflineJetPhi_,      "bTagCSVOfflineJetPhi[bTagCSVOfflineJetCounter]/F")   ;
   outTree_->Branch("bTagCSVOfflineJetEt",       bTagCSVOfflineJetEt_,       "bTagCSVOfflineJetEt[bTagCSVOfflineJetCounter]/D")    ;
   outTree_->Branch("bTagCSVOfflineJetMass",     bTagCSVOfflineJetMass_,     "bTagCSVOfflineJetMass[bTagCSVOfflineJetCounter]/D")  ;
   outTree_->Branch("bTagCSVOffline",            bTagCSVOffline_,            "bTagCSVOffline[bTagCSVOfflineJetCounter]/F")         ;    
}

bbtoDijetAnalyzer::~bbtoDijetAnalyzer()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)  
   
   // jet-related information
   // CaloJet
   delete[]  caloJetPt_                    ;   
   delete[]  caloJetEta_                   ;   
   delete[]  caloJetPhi_                   ;   
   delete[]  caloJetEt_                    ;   
   delete[]  caloJetMass_                  ;   
   delete[]  caloJetEnergy_                ;
   delete[]  caloJetPx_                    ;
   delete[]  caloJetPy_                    ;
   delete[]  caloJetPz_                    ;
   // PFJet
   delete[]  pfJetPt_                      ;       
   delete[]  pfJetEta_                     ;       
   delete[]  pfJetPhi_                     ;       
   delete[]  pfJetEt_                      ;       
   delete[]  pfJetMass_                    ;       
   delete[]  pfJetEnergy_                  ;
   delete[]  pfJetPx_                      ;
   delete[]  pfJetPy_                      ;
   delete[]  pfJetPz_                      ;
   // b-jet-tagging information
   // bTagCSVOnline       
   delete[]  bTagCSVOnlineJetPt_           ;
   delete[]  bTagCSVOnlineJetEta_          ;
   delete[]  bTagCSVOnlineJetPhi_          ;
   delete[]  bTagCSVOnlineJetEt_           ;
   delete[]  bTagCSVOnlineJetMass_         ;
   delete[]  bTagCSVOnline_                ;
   // bTagCSVOffline       
   delete[]  bTagCSVOfflineJetPt_          ;
   delete[]  bTagCSVOfflineJetEta_         ;
   delete[]  bTagCSVOfflineJetPhi_         ;
   delete[]  bTagCSVOfflineJetEt_          ;
   delete[]  bTagCSVOfflineJetMass_        ;
   delete[]  bTagCSVOffline_               ;      
}
//
// member functions
// 

// ----------- method called in the body of analyze to clear memory ----------------
// -----------      does not clear  memory for trigger results!     ----------------

void
bbtoDijetAnalyzer::clear()
{
   // set memory for branch variables
   // jet-related information
   // CaloJet 
   caloJetCounter_                     = 0                                             ;
   std::memset(caloJetPt_,            '\0',  kMaxCaloJet_          *  sizeof(double))  ;           
   std::memset(caloJetEta_,           '\0',  kMaxCaloJet_          *  sizeof(float))   ;    
   std::memset(caloJetPhi_,           '\0',  kMaxCaloJet_          *  sizeof(float))   ;   
   std::memset(caloJetEt_,            '\0',  kMaxCaloJet_          *  sizeof(double))  ;   
   std::memset(caloJetMass_,          '\0',  kMaxCaloJet_          *  sizeof(float))   ;
   std::memset(caloJetEnergy_,        '\0',  kMaxCaloJet_          *  sizeof(double))  ;
   std::memset(caloJetPx_,            '\0',  kMaxCaloJet_          *  sizeof(double))  ;
   std::memset(caloJetPy_,            '\0',  kMaxCaloJet_          *  sizeof(double))  ;
   std::memset(caloJetPz_,            '\0',  kMaxCaloJet_          *  sizeof(double))  ;
   // CaloJet   
   pfJetCounter_                       = 0                                             ;
   std::memset(pfJetPt_,              '\0',  kMaxPFJet_            *  sizeof(double))  ;           
   std::memset(pfJetEta_,             '\0',  kMaxPFJet_            *  sizeof(float))   ;    
   std::memset(pfJetPhi_,             '\0',  kMaxPFJet_            *  sizeof(float))   ;   
   std::memset(pfJetEt_,              '\0',  kMaxPFJet_            *  sizeof(double))  ;   
   std::memset(pfJetMass_,            '\0',  kMaxPFJet_            *  sizeof(double))  ;    
   std::memset(pfJetEnergy_,          '\0',  kMaxPFJet_            *  sizeof(double))  ;
   std::memset(pfJetPx_,              '\0',  kMaxPFJet_            *  sizeof(double))  ;
   std::memset(pfJetPy_,              '\0',  kMaxPFJet_            *  sizeof(double))  ;
   std::memset(pfJetPz_,              '\0',  kMaxPFJet_            *  sizeof(double))  ;
   // b-jet-tagging information
   // bTagCSVOnline
   bTagCSVOnlineJetCounter_           = 0                                              ;
   std::memset(bTagCSVOnlineJetPt_,   '\0',  kMaxBTagCSVOnline_    *  sizeof(double))  ;
   std::memset(bTagCSVOnlineJetEta_,  '\0',  kMaxBTagCSVOnline_    *  sizeof(float))   ;
   std::memset(bTagCSVOnlineJetPhi_,  '\0',  kMaxBTagCSVOnline_    *  sizeof(float))   ;
   std::memset(bTagCSVOnlineJetEt_,   '\0',  kMaxBTagCSVOnline_    *  sizeof(double))  ;
   std::memset(bTagCSVOnlineJetMass_, '\0',  kMaxBTagCSVOnline_    *  sizeof(double))  ;
   std::memset(bTagCSVOnline_,        '\0',  kMaxBTagCSVOnline_    *  sizeof(float))   ;           
   // bTagCSVOffline
   bTagCSVOfflineJetCounter_          = 0                                              ;
   std::memset(bTagCSVOfflineJetPt_,  '\0',  kMaxBTagCSVOffline_   *  sizeof(double))  ;
   std::memset(bTagCSVOfflineJetEta_, '\0',  kMaxBTagCSVOffline_   *  sizeof(float))   ;
   std::memset(bTagCSVOfflineJetPhi_, '\0',  kMaxBTagCSVOffline_   *  sizeof(float))   ;
   std::memset(bTagCSVOfflineJetEt_,  '\0',  kMaxBTagCSVOffline_   *  sizeof(double))  ;
   std::memset(bTagCSVOfflineJetMass_,'\0',  kMaxBTagCSVOffline_   *  sizeof(double))  ;
   std::memset(bTagCSVOffline_,       '\0',  kMaxBTagCSVOffline_   *  sizeof(float))   ;           
}
// ------------ method called for each event  ------------
void
bbtoDijetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{ 
    // set handles
    edm::Handle<edm::TriggerResults>       triggerResultsHandle    ;
    edm::Handle<trigger::TriggerEvent>     triggerEventHandle      ;
    edm::Handle<reco::CaloJetCollection>   caloJetsHandle          ;
    edm::Handle<reco::PFJetCollection>     pfJetsHandle            ;
    edm::Handle<reco::JetTagCollection>    bTagCSVOnlineHandle     ;
    edm::Handle<reco::JetTagCollection>    bTagCSVOfflineHandle    ; 
  
    // get event header
    eventHeader_.analyze(iEvent, outTree_);

    // get trigger results
    iEvent.getByToken(triggerResultsToken_, triggerResultsHandle)                        ;
    edm::TriggerNames const& triggerNames = iEvent.triggerNames(*triggerResultsHandle)   ;
    for(unsigned int itrig = 0; itrig != triggerResultsHandle->size(); ++itrig) {
      TString triggerName = triggerNames.triggerName(itrig);
      if(firstEvent_){
        outTree_->Branch(triggerName, &triggerPass_[itrig], triggerName+"/i");
     }
     bool accept = triggerResultsHandle->accept(itrig);
     if(accept) triggerPass_[itrig] = 1;
     else       triggerPass_[itrig] = 0; 
    }

    if(firstEvent_) firstEvent_ = false; 

    // if the required collections are available, fill the corresponding tree branches
    // get jet-related information
    // CaloJets
    if(iEvent.getByToken(caloJetsToken_, caloJetsHandle)){
      const reco::CaloJetCollection & caloJets = *(caloJetsHandle.product());
      caloJetCounter_ = caloJets.size();
      for(int i = 0; i != caloJetCounter_; ++i) {
        caloJetPt_[i]     = caloJets[i].pt()   ;
        caloJetEta_[i]    = caloJets[i].eta()  ;
        caloJetPhi_[i]    = caloJets[i].phi()  ;
        caloJetEt_[i]     = caloJets[i].et()   ;
        caloJetMass_[i]   = caloJets[i].mass() ;
        caloJetEnergy_[i] = caloJets[i].energy()  ;
        caloJetPx_[i]     = caloJets[i].px()   ;
        caloJetPy_[i]     = caloJets[i].py()   ;
        caloJetPz_[i]     = caloJets[i].pz()   ;
      }
    }
    // PFJets
    if(iEvent.getByToken(pfJetsToken_, pfJetsHandle)){
      const reco::PFJetCollection & pfJets = *(pfJetsHandle.product());
      pfJetCounter_ = pfJets.size();
      for(int i = 0; i != pfJetCounter_; ++i) {
        pfJetPt_[i]     = pfJets[i].pt()   ;
        pfJetEta_[i]    = pfJets[i].eta()  ;
        pfJetPhi_[i]    = pfJets[i].phi()  ;
        pfJetEt_[i]     = pfJets[i].et()   ;
        pfJetMass_[i]   = pfJets[i].mass() ;
        pfJetEnergy_[i] = pfJets[i].energy() ;
        pfJetPx_[i]     = pfJets[i].px()   ;
        pfJetPy_[i]     = pfJets[i].py()   ;
        pfJetPz_[i]     = pfJets[i].pz()   ;
      }
    }
    // get b-jet-tagging information
    // bTagCSVOnline
    if(iEvent.getByToken(bTagCSVOnlineToken_, bTagCSVOnlineHandle)){
      const reco::JetTagCollection & bTagCSVOnline = *(bTagCSVOnlineHandle.product());
      bTagCSVOnlineJetCounter_ = bTagCSVOnline.size();
      for(int i = 0; i != bTagCSVOnlineJetCounter_; ++i){
        // save the tag and corresponding pt
        // see  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagEdAnalyzer17X
        bTagCSVOnlineJetPt_[i]   = bTagCSVOnline[i].first->pt()    ;
        bTagCSVOnlineJetEta_[i]  = bTagCSVOnline[i].first->eta()   ;
        bTagCSVOnlineJetPhi_[i]  = bTagCSVOnline[i].first->phi()   ;
        bTagCSVOnlineJetEt_[i]   = bTagCSVOnline[i].first->et()    ;
        bTagCSVOnlineJetMass_[i] = bTagCSVOnline[i].first->mass()  ;
        bTagCSVOnline_[i]        = bTagCSVOnline[i].second         ;
      }
    }
    // bTagCSVOffline         
    if(iEvent.getByToken(bTagCSVOfflineToken_, bTagCSVOfflineHandle)){
      const reco::JetTagCollection & bTagCSVOffline = *(bTagCSVOfflineHandle.product());
      bTagCSVOfflineJetCounter_ = bTagCSVOffline.size();
      for(int i = 0; i != bTagCSVOfflineJetCounter_; ++i){
        // save the tag and corresponding pt
        // see  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagEdAnalyzer17X
        bTagCSVOfflineJetPt_[i]   = bTagCSVOffline[i].first->pt()  ;
        bTagCSVOfflineJetEta_[i]  = bTagCSVOffline[i].first->eta() ;
        bTagCSVOfflineJetPhi_[i]  = bTagCSVOffline[i].first->phi() ;
        bTagCSVOfflineJetEt_[i]   = bTagCSVOffline[i].first->et()  ;
        bTagCSVOfflineJetMass_[i] = bTagCSVOffline[i].first->mass();
        bTagCSVOffline_[i]        = bTagCSVOffline[i].second       ;
      }
    }

    outTree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
bbtoDijetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbtoDijetAnalyzer::endJob() 
{
}
// ------------ method called when starting to processes a run  ------------
void 
bbtoDijetAnalyzer::beginRun(edm::Run const &run, edm::EventSetup const &es)
{
    bool changed;

    if (!hltConfig_.init(run, es, processName_, changed)) {
      edm::LogError("bbtoDijetAnalyzer") << "Initialization of HLTConfigProvider failed.";
      return;
    }
}
// ------------ method called when starting to processes a run  ------------
void
bbtoDijetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bbtoDijetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("processName", "TEST");
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults::TEST"));
  desc.add<edm::InputTag>("triggerEvent",   edm::InputTag("hltTriggerSummaryAOD", "", "TEST")); 
  desc.add<edm::InputTag>("caloJets",       edm::InputTag("ak4CaloJets::RECO"));
  desc.add<edm::InputTag>("pfJets",         edm::InputTag("ak4PFJets::RECO"));
  desc.add<edm::InputTag>("bTagCSVOnline",  edm::InputTag("hltCombinedSecondaryVertexBJetTagsCalo"));
  desc.add<edm::InputTag>("bTagCSVOffline", edm::InputTag("pfCombinedSecondaryVertexV2BJetTags"));
  descriptions.add("bbtoDijetAnalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbtoDijetAnalyzer);
