// -*- C++ -*-
//
// Package:    HLTrigger/SimpleHLTAnalyzer
// Class:      SimpleHLTAnalyzer
// 
/**\class SimpleHLTAnalyzer SimpleHLTAnalyzer.cc HLTrigger/SimpleHLTAnalyzer/plugins/SimpleHLTAnalyzer.cc

 Description: save trigger decisions into a TTree, as well as physics objects for computing efficiencies

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Valentina Dutta
//         Created:  Tue, 03 Mar 2015 10:42:08 GMT
//
// Edited by: 	     Arthur Avkhadiev
// Began modifying:  Wed Jul  6 14:03:58 CEST 2016


#include <memory>

#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h" 
#include "DataFormats/Luminosity/interface/LumiDetails.h" 
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class SimpleHLTAnalyzer : public edm::EDAnalyzer {
  public:
    explicit SimpleHLTAnalyzer(const edm::ParameterSet&); 
    ~SimpleHLTAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

    HLTConfigProvider hltConfig_;
    std::string hltProcess_;

    // Tokens for trigger results, kinematics, and btags
    edm::EDGetTokenT<edm::TriggerResults>    trigresultsToken_;
    edm::EDGetTokenT<trigger::TriggerEvent>  trigsummaryToken_;
    edm::EDGetTokenT<reco::PFJetCollection>  pfjetsToken_; 
    edm::EDGetTokenT<reco::JetTagCollection> btagsCSVOnlineToken_; 
    edm::EDGetTokenT<reco::JetTagCollection> btagsCSVOfflineToken_;
    edm::EDGetTokenT<LumiSummary>            lumiToken_;
    
    unsigned int run_;
    unsigned int lumi_;
    unsigned int evt_;
    double AvgInstDelLumi_;
 
    bool firstEvent_;

    const int maxResults_;

    unsigned int *passtrig_;
    // b-tagging discriminants
    int nbtagCSVOnlinejets_;
    int nbtagCSVOfflinejets_;
    float *bTagCSVOnline_;
    float *bTagCSVOffline_;
    // PFJet objects
    int npfjets_;
    float *pfjetpt_, *pfjeteta_, *pfjetphi_, *pfjetmass_;

    TTree* outTree_;

};

SimpleHLTAnalyzer::SimpleHLTAnalyzer(const edm::ParameterSet& iConfig) :
  // the assigned strings have to match the names in the config file  
  hltProcess_          (iConfig.getParameter<std::string>("hltProcess")),
  trigresultsToken_    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigResults"))),
  trigsummaryToken_    (consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trigSummary"))),
  pfjetsToken_         (consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJets"))),
  btagsCSVOnlineToken_ (consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("bTagsCSVOnline"))),
  btagsCSVOfflineToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("bTagsCSVOffline"))), 
  lumiToken_           (consumes<LumiSummary,edm::InLumi>(iConfig.getParameter<edm::InputTag>("lumiProducer"))),
  run_  (0),
  lumi_ (0),
  evt_  (0),
  AvgInstDelLumi_ (-999.), 
  firstEvent_ (true),
  maxResults_ (1000)
{

  passtrig_ = new unsigned int[maxResults_];
  AvgInstDelLumi_ = -999.;
  bTagCSVOnline_ = new float[maxResults_];
  bTagCSVOffline_ = new float[maxResults_];
  nbtagCSVOnlinejets_ = 0;
  nbtagCSVOfflinejets_ = 0;
  npfjets_ = 0;
  pfjetpt_ = new float[maxResults_];
  pfjeteta_ = new float[maxResults_];
  pfjetphi_ = new float[maxResults_];
  pfjetmass_ = new float[maxResults_];

  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("HLTAnalysis","");

// specifying the branch: slash letter is ROOT's way for figuring out types:
// I is int
// F is float
// D is double
  outTree_->Branch("run",  &run_,  "run/i");
  outTree_->Branch("lumi", &lumi_, "lumi/i");
  outTree_->Branch("AvgInstDelLumi", &AvgInstDelLumi_, "AvgInstDelLumi/D");
  outTree_->Branch("evt",  &evt_,  "evt/i");
  outTree_->Branch("nbtagCSVOnlinejets", &nbtagCSVOnlinejets_, "nbtagCSVOnlinejets/I");
  outTree_->Branch("bTagCSVOnline",  bTagCSVOnline_, "bTagCSVOnline[nbtagCSVOnlinejets]/F"); 
  outTree_->Branch("nbtagCSVOfflinejets", &nbtagCSVOfflinejets_, "nbtagCSVOfflinejets/I");
  outTree_->Branch("bTagCSVOffline",  bTagCSVOffline_, "bTagCSVOffline[nbtagCSVOfflinejets]/F");
  outTree_->Branch("npfjets",  &npfjets_,  "npfjets/I");
  outTree_->Branch("pfjetpt",  pfjetpt_,   "pfjetpt[npfjets]/F");
  outTree_->Branch("pfjeteta", pfjeteta_,   "pfjeteta[npfjets]/F");
  outTree_->Branch("pfjetphi", pfjetphi_,  "pfjetphi[npfjets]/F");
  outTree_->Branch("pfjetmass", pfjetmass_, "pfjetmass[npfjets]/F");
}

SimpleHLTAnalyzer::~SimpleHLTAnalyzer()
{}

// ------------ method called for each event  ------------
void
SimpleHLTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // event header
  run_  = iEvent.id().run();
  lumi_ = iEvent.luminosityBlock();
  evt_  = iEvent.id().event();

  // average ins delivered luminosity
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_1_4/doc/html/d9/d79/EventHeader_8cc_source.html
  bool lumiException = false;
  const edm::LuminosityBlock& iLumi = iEvent.getLuminosityBlock();
  edm::Handle<LumiSummary> lumiSummary; 
  try{
    iLumi.getByToken(lumiToken_, lumiSummary);
    lumiSummary->isValid();
  }
  catch(cms::Exception&){
    lumiException = true;
  }
  if(!lumiException)
      AvgInstDelLumi_ = lumiSummary->avgInsDelLumi();
  else
      AvgInstDelLumi_ = -999.;

  // trigger results
  edm::Handle<edm::TriggerResults> trigresults;
  iEvent.getByToken(trigresultsToken_, trigresults);

  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*trigresults); 

  for(unsigned int itrig = 0; itrig < trigresults->size(); ++itrig) {
    TString trigname = triggerNames.triggerName(itrig);
    if(firstEvent_) {
      outTree_->Branch(trigname,&passtrig_[itrig],trigname+"/i");
    }
    bool accept = trigresults->accept(itrig);
    if(accept) passtrig_[itrig] = 1;
    else       passtrig_[itrig] = 0;
  }

  if(firstEvent_) firstEvent_ = false;

  // kinematics from PFJets
  edm::Handle<reco::PFJetCollection> pfjetHandle;
  iEvent.getByToken(pfjetsToken_, pfjetHandle);
  const reco::PFJetCollection & pfjets = *(pfjetHandle.product());

  npfjets_ = pfjets.size();
  
  // loop over jets
  for(int i = 0; i != npfjets_; ++i) {
    pfjetpt_[i]   = pfjets[i].pt();
    pfjeteta_[i]  = pfjets[i].eta();
    pfjetphi_[i]  = pfjets[i].phi();
    pfjetmass_[i] = pfjets[i].mass();
  }

  // b-tagging discriminants
  
  // CSVOnline
  // see https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/de/d43/HLTScoutingCaloProducer_8cc_source.html
  edm::Handle<reco::JetTagCollection> btagCSVOnlineHandle;
  bool haveJetBTags = false;
  if(iEvent.getByToken(btagsCSVOnlineToken_, btagCSVOnlineHandle)){
    haveJetBTags = true; 
   }

  if(haveJetBTags){
    const reco::JetTagCollection & btagsCSVOnline = *(btagCSVOnlineHandle.product());
    nbtagCSVOnlinejets_ = btagsCSVOnline.size();
    //  loop over jets 
    for(int i = 0; i != nbtagCSVOnlinejets_; ++i) {
        bTagCSVOnline_[i] = btagsCSVOnline[i].second; // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagEdAnalyzer17X
    }
  }
  
  // CSVOffline
  edm::Handle<reco::JetTagCollection> btagCSVOfflineHandle;
  iEvent.getByToken(btagsCSVOfflineToken_, btagCSVOfflineHandle);
  const reco::JetTagCollection & btagsCSVOffline = *(btagCSVOfflineHandle.product());

  nbtagCSVOfflinejets_ = btagsCSVOffline.size();
  //  loop over jets 
  for(int i = 0; i != nbtagCSVOfflinejets_; ++i) {
    bTagCSVOffline_[i] = btagsCSVOffline[i].second; // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookBTagEdAnalyzer17X
  }

  outTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
SimpleHLTAnalyzer::beginJob()
{}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleHLTAnalyzer::endJob() 
{}

// ------------ method called when starting to processes a run  ------------
void 
SimpleHLTAnalyzer::beginRun(edm::Run const &run, edm::EventSetup const &es)
{

  bool changed;

  if (!hltConfig_.init(run, es, hltProcess_, changed)) {
    edm::LogError("SimpleHLTAnalyzer") << "Initialization of HLTConfigProvider failed!!";
    return;
  }

}

// ------------ method called when ending the processing of a run  ------------
void 
SimpleHLTAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleHLTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleHLTAnalyzer);
