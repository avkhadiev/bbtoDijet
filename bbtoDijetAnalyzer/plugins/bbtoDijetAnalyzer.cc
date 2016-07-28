// -*- C++ -*-
//
// Package:    bbtoDijet/bbtoDijetAnalyzer
// Class:      bbtoDijetAnalyzer
// 
/**\class bbtoDijetAnalyzer bbtoDijetAnalyzer.cc bbtoDijet/bbtoDijetAnalyzer/plugins/bbtoDijetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Artur Avkhadiev
//         Created:  Fri, 22 Jul 2016 08:05:25 GMT
//
//


// system include files
// uncluded by default
#include <memory>
// added in
#include <TTree.h>

// user include files
// included by default
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// added in
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/dd/d5c/HLTBitAnalyzer_8h_source.html
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/d04/HLTJets_8h_source.html
#include "FWCore/Framework/interface/ConsumesCollector.h"
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/df/d03/HLTEventAnalyzerRAW_8cc_source.html
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d5/d74/HLTBitAnalyzer_8cc_source.html
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/dd/d5c/HLTBitAnalyzer_8h_source.html
#include "HLTrigger/HLTanalyzers/interface/EventHeader.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/dc0/HLTEventAnalyzerAOD_8cc_source.html
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/Handle.h"

// event content access
// event-header information
#include "DataFormats/Luminosity/interface/LumiSummary.h" 
#include "DataFormats/Luminosity/interface/LumiDetails.h"

// trigger-related information
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h" 
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

// jet-related information
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/d04/HLTJets_8h_source.html
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// b-jet-tagging information
// from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_13/doc/html/d2/d2e/HLTBJet_8h_source.html
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

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
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void clear()  override;
      virtual void endJob() override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

// input variables
// analysis tools
HLTConfigProvider  hltConfig_    ;
EventHeader        eventHeader_  ;
std::string        processName_  ;

// input tags
const edm::InputTag                                     triggerResultsTag_            ;
const edm::InputTag                                     triggerEventTag_              ;

const edm::InputTag                                     pfJetsTag_                    ;
const edm::InputTag                                     caloJetsTag_                  ;
const edm::InputTag                                     bTagCSVOnlineTag_             ; 
const edm::InputTag                                     bTagCSVOfflineTag_            ;

// tokens                                              
const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken_          ;
const edm::EDGetTokenT<trigger::TriggerEvent>           triggerEventToken_            ;
const edm::EDGetTokemT<reco::CaloJetCollecton>          caloJetsToken_                ;
const edm::EDGetTokenT<reco::PFJetCollection>           pfJetsToken_                  ;
const edm::EDGetTokenT<reco::JetTagCollection>          bTagCSVOnlineToken_           ; 
const edm::EDGetTokenT<reco::JetTagCollection>          bTagCSVOfflineToken_          ;

// tree variables
// array sizes
static const int kMaxTriggerPass_    ;
static const int kMaxCaloJet_        ;
static const int kMaxPFJet_          ;
static const int kMaxBTagCSVOnline_  ;
static const int kMaxBTagCSVOffline_ ;

// set of variables for 
// event-header information: included in eventHeader_
// trigger-related information
bool          firstEvent_                 ;   // will create new branch for trigger path if firstEvent_ = true
unsigned int  triggerPass_                ;
// jet-related information
// CaloJet
int           caloJetCounter_             ;     
float        *caloJetPt_                  ;  
float        *caloJetEta_                 ; 
float        *caloJetPhi_                 ;  
float        *caloJetEt_                  ;  
float        *caloJetMass_                ; 
// PFJet
int           pfJetCounter_               ;     
float        *pfJetPt_                    ;  
float        *pfJetEta_                   ; 
float        *pfJetPhi_                   ;  
float        *pfJetEt_                    ;  
float        *pfJetMass_                  ; 
// b-jet-tagging information
// bTagCSVOnline
int           bTagCSVOnlineJetCounter_    ;     
float        *bTagCSVOnlineJetPt_         ;
float        *bTagCSVOnline_              ;
// bTagCSVOffline
int           bTagCSVOfflineJetCounter_   ;
float        *bTagCSVOfflineJetPt_        ;
float        *bTagCSVOffline_             ;

TTree* outTree_;

// constructors and destructor
//
bbtoDijetAnalyzer::bbtoDijetAnalyzer(const edm::ParameterSet& iConfig) :
  processName_                       (iConfig.getParameter<std::string>("processName"), 
  eventHeader_                       (),                                              ,       
  // initialize input tags 
  // (strings should match variable names in the *.py config file)
  triggerResultsTag_                 (iConfig.getParameter<std::string>)("triggerResults"), 
  triggerEventTag_                   (iConfig.getParameter<std::string>)("triggerEvent"),                     
  pfJetsTag_                         (iConfig.getParameter<std::string>)("caloJets"),            
  caloJetsTag_                       (iConfig.getParameter<std::string>)("pfJets"),            
  bTagCSVOnlineTag_                  (iConfig.getParameter<std::string>)("bTagCSVOnline"),            
  bTagCSVOfflineTag_                 (iConfig.getParameter<std::string>)("bTagCSVOffline"),            
  // initialize tokens
  triggerResultsToken_               (consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerEventToken_                 (consumes<trigger::TriggerEvent>(triggerEventTag_)), 
  caloJetsToken_                     (consumes<reco::CaloJetCollection>(caloJetsTag_)),
  pfJetsToken_                       (consumes<reco::PFJetCollection>(pfJetsTag_)),
  bTagCSVOnlineToken_                (consumes<reco::JetTagCollection>(bTagCSVOnlineTag_)),
  bTagCSVOfflineToken_               (consumes<reco::JetTagCollection>(bTagCSVOfflineTag_)), 
  // initialize tree variables
  firstEvent_                        (true),
  kMaxTriggerPass_                   (10000),
  kMaxCaloJet_                       (10000),
  kMaxPFJet_                         (10000),
  kMaxBTagCSVOnline_                 (1000),
  kMaxBTagCSVOffline_                (1000)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
  
   /*  Setup the analysis to put the branch-variables into the tree. */

   // set of variables for
   // event-header information: included in eventHeader_ 
   // trigger-related information
   triggerPass_               = new unsigned int[kMaxTriggerPass]   ;
   // jet-related information
   // CaloJet
   caloJetCounter_            = 0                                   ;
   caloJetPt_                 = new float[kMaxCaloJet]              ;        
   caloJetEta_                = new float[kMaxCaloJet]              ; 
   caloJetPhi_                = new float[kMaxCaloJet]              ; 
   caloJetEt_                 = new float[kMaxCaloJet]              ; 
   caloJetMass_               = new float[kMaxCaloJet]              ; 
   // PFJet
   pfJetCounter_              = 0                                   ;
   pfJetPt_                   = new float[kMaxPFJet]                ;         
   pfJetEta_                  = new float[kMaxPFJet]                ;  
   pfJetPhi_                  = new float[kMaxPFJet]                ;  
   pfJetEt_                   = new float[kMaxPFJet]                ;  
   pfJetMass_                 = new float[kMaxPFJet]                ; 
   // b-jet-tagging information
   // bTagCSVOnline
   bTagCSVOnlineJetCounter_   = 0                                   ;
   bTagCSVOnlineJetPt_        = new float[kMaxBTagCSVOnline]        ;
   bTagCSVOnline_             = new float[kMaxBTagCSVOnline]        ;
   // bTagCSVOffline 
   bTagCSVOfflineJetCounter_  = 0                                   ;  
   bTagCSVOfflineJetPt_       = new float[kMaxBTagCSVOffline]       ;
   bTagCSVOffline_            = new float[kMaxBTagCSVOffline]       ;

   // open the tree file and initialize the tree
   edm::Service<TFileService> fs                                    ;
   outTree_ = fs->make<TTree>("efficiencyTree", "")                 ;

   // setup event header analysis
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
   outTree_->Branch("caloJetCounter",           &caloJetCounter_,            "caloJetCounter/I")                            ;
   outTree_->Branch("caloJetPt",                 caloJetPt_,                 "caloJetPT[caloJetCounter]/F")                 ;    
   outTree_->Branch("caloJetEta",                caloJetEta_,                "caloJetEta[caloJetCounter]/F")                ;    
   outTree_->Branch("caloJetPhi",                caloJetPhi_,                "caloJetPhi[caloJetCounter]/F")                ;    
   outTree_->Branch("caloJetEt",                 caloJetEt_,                 "caloJetEt[caloJetCounter]/F")                 ;    
   outTree_->Branch("caloJetMass",               caloJetMass_,               "caloJetMass[caloJetCounter]/F")               ;    
   // PFJet
   outTree_->Branch("pfJetCounter",             &pfJetCounter_,              "pfJetCounter/I")                              ;
   outTree_->Branch("pfJetPt",                   pfJetPt_,                   "pfJetPT[pfJetCounter]/F")                     ;    
   outTree_->Branch("pfJetEta",                  pfJetEta_,                  "pfJetEta[pfJetCounter]/F")                    ;    
   outTree_->Branch("pfJetPhi",                  pfJetPhi_,                  "pfJetPhi[pfJetCounter]/F")                    ;    
   outTree_->Branch("pfJetEt",                   pfJetEt_,                   "pfJetEt[pfJetCounter]/F")                     ;    
   outTree_->Branch("pfJetMass",                 pfJetMass_,                 "pfJetMass[pfJetCounter]/F")                   ;    
   // b-jet-tagging information
   // bTagCSVOnline
   outTree_->Branch("bTagCSVOnlineJetCounter",  &bTagCSVOnlineJetCounter_,   "bTagCSVOnlineJetCounter/I")                   ;
   outTree_->Branch("bTagCSVOnlineJetPt",        bTagCSVOnlineJetPt_,        "bTagCSVOnlineJetPt[bTagCSVOnlineJetCounter]/F")                        ;
   outTree_->Branch("bTagCSVOnline",             bTagCSVOnline_,             "bTagCSVOnline[bTagCSVOnlineJetCounter]/F")    ;    
   // bTagCSVOffline
   outTree_->Branch("bTagCSVOfflineJetCounter", &bTagCSVOfflineJetCounter_,  "bTagCSVOfflineJetCounter/I")                  ;
   outTree_->Branch("bTagCSVOfflineJetPt",       bTagCSVOfflineJetPt_,       "bTagCSVOfflineJetPt[bTagCSVOfflineJetCounter]/F")                       ;
   outTree_->Branch("bTagCSVOffline",            bTagCSVOffline_,            "bTagCSVOffline[bTagCSVOfflineJetCounter]/F")  ;    
}

bbtoDijetAnalyzer::~bbtoDijetAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)  
   
   // jet-related information
   // CaloJet
   delete[]  caloJetCounter_               ;
   delete[]  caloJetPt_                    ;   
   delete[]  caloJetEta_                   ;   
   delete[]  caloJetPhi_                   ;   
   delete[]  caloJetEt_                    ;   
   delete[]  caloJetMass_                  ;   
   // PFJet
   delete[]  pfJetCounter_                 ;       
   delete[]  pfJetPt_                      ;       
   delete[]  pfJetEta_                     ;       
   delete[]  pfJetPhi_                     ;       
   delete[]  pfJetEt_                      ;       
   delete[]  pfJetMass_                    ;       
   // b-jet-tagging information
   // bTagCSVOnline       
   delete[]  bTagCSVOnlineJetCounter_      ;
   delete[]  bTagCSVOnlineJetPt_           ;
   delete[]  bTagCSVOnline_                ;
   // bTagCSVOffline       
   delete[]  bTagCSVOfflineJetCounter_     ;
   delete[]  bTagCSVOfflineJetPt_          ;
   delete[]  bTagCSVOffline_               ;      
   
//
// member functions
// 

// ------------ method called for each event  ------------
void
bbtoDijetAnalyzer::clear()
{
   // set memory for branch variables
   // jet-related information
   // CaloJet 
   caloJetCounter_                    = 0                                            ;
   std::memset(caloJetPt_,           '\0',  kMaxCaloJet_          *  sizeof(float))  ;           
   std::memset(caloJetEta_,          '\0',  kMaxCaloJet_          *  sizeof(float))  ;    
   std::memset(caloJetPhi_,          '\0',  kMaxCaloJet_          *  sizeof(float))  ;   
   std::memset(caloJetEt_,           '\0',  kMaxCaloJet_          *  sizeof(float))  ;   
   std::memset(caloJetMass_,         '\0',  kMaxCaloJet_          *  sizeof(float))  ;    
   // CaloJet   
   pfJetCounter_                      = 0                                            ;
   std::memset(pfJetPt_,             '\0',  kMaxPFJet_            *  sizeof(float))  ;           
   std::memset(pfJetEta_,            '\0',  kMaxPFJet_            *  sizeof(float))  ;    
   std::memset(pfJetPhi_,            '\0',  kMaxPFJet_            *  sizeof(float))  ;   
   std::memset(pfJetEt_,             '\0',  kMaxPFJet_            *  sizeof(float))  ;   
   std::memset(pfJetMass_,           '\0',  kMaxPFJet_            *  sizeof(float))  ;    
   // b-jet-tagging information
   // bTagCSVOnline
   bTagCSVOnlineJetCounter_           = 0                                           ;
   std::memset(bTagCSVOnlineJetPt_,   '\0',  kMaxBTagCSVOffline_  *  sizeof(float)) ;
   std::memset(bTagCSVOnline_,        '\0',  kMaxBTagCSVOffline_  *  sizeof(float)) ;           
   // bTagCSVOffline
   bTagCSVOfflineJetCounter_          = 0                                            ;
   std::memset(bTagCSVOfflineJetPt_, '\0',  kMaxBTagCSVOffline_   *  sizeof(float))  ;
   std::memset(bTagCSVOffline_,      '\0',  kMaxBTagCSVOffline_   *  sizeof(float))  ;           
}

void
bbtoDijetAnalyzer::analyzeJets(const edm::View<reco::Jet>)

void
bbtoDijetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{ 
    // set handles
    const edm::Handle<edm::TriggerResults>       triggerResultsHandle    ;
    const edm::Handle<trigger:TriggerEvent>      triggerEventHandle      ;
    const edm::Handle<reco::CaloJetCollection>   caloJetsHandle          ;
    const edm::Handle<reco::PFJetCollection>     pfJetsHandle            ;
    const edm::Handle<reco::JetTagCollection>    bTagCSVOnlineHandle     ;
    const edm::Handle<reco::JetTagCollection>    bTagCSVOfflineHandle    ; 
  
    // get trigger results
    iEvent.getByToken(triggerResultsToken_, triggerResultsHandle)                        ;
    edm::TriggerNames const& triggerNames = iEvent.triggerNames( *triggerResultsHandle ) ;
    for(unsigned int itrig = ; triggerResultsHandle->size(); ++itrig) {
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
    if(iEvent.getByToken(caloJetsTokem_, caloJetsHandle)){
      const reco::CaloJetCollection & caloJets = *(caloJetsHandle.product());
      caloJetCounter_ = caloJets.size();
      for(int i = 0; i != caloJetCounter_; ++i) {
        caloJetPt_[i]   = caloJets[i].pt()   ;
        caloJetEta_[i]  = caloJets[i].eta()  ;
        caloJetPhi_[i]  = caloJets[i].phi()  ;
        caloJetEt_[i]   = caloJets[i].et()   ;
        caloJetMass_[i] = caloJets[i].mass() ;
      }
    }
    // PFJets
    if(iEvent.getByToken(pfJetsTokem_, pfJetsHandle)){
      const reco::PFJetCollection & pfJets = *(pfJetsHandle.product());
      pfJetCounter_ = pfJets.size();
      for(int i = 0; i != pfJetCounter_; ++i) {
        pfJetPt_[i]     = pfJets[i].pt()   ;
        pfJetEta_[i]    = pfJets[i].eta()  ;
        pfJetPhi_[i]    = pfJets[i].phi()  ;
        pfJetEt_[i]     = pfJets[i].et()   ;
        pfJetMass_[i]   = pfJets[i].mass() ;
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
        bTagCSVOnlineJetPt_[i] = bTagCSVOnline[i].first;
        bTagCSVOnline_[i]      = bTagCSVOnline[i].second;
      }
    }
    // bTagCSVOffline         
    if(iEvent.getByToken(bTagCSVOfflineToken_, bTagCSVOfflineHandle)){
      const reco::JetTagCollection & bTagCSVOffline = *(bTagCSVOfflineHandle.product());
      bTagCSVOfflineJetCounter_ = bTagCSVOffline.size();
      for(int i = 0; i != bTagCSVOfflineJetCounter_; ++i){
        // save the tag and corresponding pt
        // see  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagEdAnalyzer17X
        bTagCSVOfflineJetPt_[i] = bTagCSVOffline[i].first;
        bTagCSVOffline_[i]      = bTagCSVOffline[i].second;
      }
    }

    outTree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
DijetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bbtoDijetAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bbtoDijetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("processName", "TEST");
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults::TEST"));
  desc.add<edm::InputTag>("triggerEvent",   edm::InputTag("hltTriggerSummaryAOD", "", "TEST");
  desc.add<edm::InputTag>("pfJets",         edm::InputTag(""));
  desc.add<edm::InputTag>("caloJets",       edm::InputTag(""));
  desc.add<edm::InputTag>("bTagCSVOnline",  edm::InputTag("hltCombinedSecondaryVertexBJetTagsCalo"));
  desc.add<edm::InputTag>("bTagCSVOffline", edm::InputTag("pfCombinedSecondaryVertexV2BJetTags"));
  descriptions.add("bbtoDijetAnalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bbtoDijetAnalyzer);
