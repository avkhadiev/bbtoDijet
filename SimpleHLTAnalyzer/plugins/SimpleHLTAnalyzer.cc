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
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"


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

    // Filter for accessing online kinematics
    // edm::InputTag trigFilterJet_;

    // Tokens for trigger results, kinematics, and btags
    edm::EDGetTokenT<edm::TriggerResults>    trigresultsToken_;
    edm::EDGetTokenT<trigger::TriggerEvent>  trigsummaryToken_;
    edm::EDGetTokenT<reco::PFJetCollection> pfjetsToken_; 
    edm::EDGetTokenT<reco::JetTagCollection> btagsCSVOnlineToken_; 
    edm::EDGetTokenT<reco::JetTagCollection> btagsCSVOfflineToken_;


    
// tree variablesunsigned int run_;
    unsigned int lumi_;
    unsigned int evt_;
 
    bool firstEvent_;

    // hard coded EDM relic:
    // create large arrays that will definitely fit 
    // all variables for a given object 
    // (e.g., pt's of a jet)
    const int maxResults_;

    unsigned int *passtrig_;
    // b-tagging discriminants
    int nbtagCSVOnlinejets_;
    int nbtagCSVOfflinejets_;
    // leading jets
    // float *bTagCSVOnline1_, *bTagCSVOnline2_, *bTagCSVOnline3_, *bTagCSVOnline4_;
    // float *bTagCSVOffline1_, *bTagCSVOffline2_, *bTagCSVOffline3_, *bTagCSVOffline4_;
    // arrays for sorting  
    float *bTagCSVOnline_;
    float *bTagCSVOffline_;
    // PFJet objects
    int npfjets_;
    // leading jets
    // float *pfjetpt1_, *pfjeteta1_, *pfjetphi1_, *pfjetmass1_;
    // float *pfjetpt2_, *pfjeteta2_, *pfjetphi2_, *pfjetmass2_;
    // float *pfjetpt3_, *pfjeteta3_, *pfjetphi3_, *pfjetmass3_;
    // float *pfjetpt4_, *pfjeteta4_, *pfjetphi4_, *pfjetmass4_;
    // arrays for sorting 
    float *pfjetpt_, *pfjeteta_, *pfjetphi_, *pfjetmass_;
    // for debugging
    int isptOrdered_;
    // online kinematics
    // int nhltjets_;
    // float *hltjetpt_, *hltjeteta_, *hltjetphi_, *hltjetmass_;

    TTree* outTree_;

};

SimpleHLTAnalyzer::SimpleHLTAnalyzer(const edm::ParameterSet& iConfig) :
  // the assigned strings have to match the names in the config file  
  hltProcess_       (iConfig.getParameter<std::string>("hltProcess")),
  // trigFilterJet_    (iConfig.getParameter<edm::InputTag>("trigFilterJet")),
  trigresultsToken_ (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigResults"))),
  trigsummaryToken_ (consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trigSummary"))),
  pfjetsToken_     (consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJets"))),
  btagsCSVOnlineToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("bTagsCSVOnline"))),
  btagsCSVOfflineToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("bTagsCSVOffline"))), 
  run_  (0),
  lumi_ (0),
  evt_  (0),
  firstEvent_ (true),
  maxResults_ (1000)
{

  passtrig_ = new unsigned int[maxResults_];
  //btagging discriminants 
  // 1st leading jet
  // bTagCSVOnline1_ = 0;
  // bTagCSVOffline1_ = 0;
  // 2nd leading jet
  // bTagCSVOnline2_ = 0;
  // bTagCSVOffline3_ = 0;
  // 3rd leading jet
  // bTagCSVOnline3_ = 0;
  // bTagCSVOffline3_ = 0;
  // 4th leading jet
  // bTagCSVOnline4_ = 0;
  // bTagCSVOffline4_ = 0;
  // btagging arrays 
  bTagCSVOnline_ = new float[maxResults_];
  bTagCSVOffline_ = new float[maxResults_];
  nbtagCSVOnlinejets_ = 0;
  nbtagCSVOfflinejets_ = 0;
  // pfjet kinematics
  npfjets_ = 0;
  // 1st leading jet
  // pfjetpt1_ = 0;
  // pfjeteta1_ = 0;
  // pfjetphi1_ = 0;
  // pfjetmass1_ = 0;
  // 2nd leading jet
  // pfjetpt2_ = 0;
  // pfjeteta2_ = 0;
  // pfjetphi2_ = 0;
  // pfjetmass2_ = 0;
  // 3rd leading jet
  // pfjetpt3_ = 0;
  //pfjeteta3_ = 0;
  // pfjetphi3_ = 0;
  // pfjetmass3_ = 0;
  // 4th leading jet
  // pfjetpt4_ = 0;
  // pfjeteta4_ = 0;
  // pfjetphi4_ = 0;
  // pfjetmass4_ = 0;
  // jet arrrays
  pfjetpt_ = new float[maxResults_];
  pfjeteta_ = new float[maxResults_];
  pfjetphi_ = new float[maxResults_];
  pfjetmass_ = new float[maxResults_];
 // for debugging 
  isptOrdered_ = 1;
 // online info
 // nhltjets_ = 0;
 // hltjetpt_ = new float[maxResults_];
 // hltjeteta_ = new float[maxResults_];
 // hltjetphi_ = new float[maxResults_];
 // hltjetmass_ = new float[maxResults_];

  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("HLTAnalysis","");

// specifying the branch: slash letter is ROOT's way for figuring out types:
// I is int
// F is float
// D is double
  outTree_->Branch("run",  &run_,  "run/i");
  outTree_->Branch("lumi", &lumi_, "lumi/i");
  outTree_->Branch("evt",  &evt_,  "evt/i");
//  outTree_->Branch("bTagCSVOnline1", bTagCSVOnline1_, "bTagCSVOnline1/F");
//  outTree_->Branch("bTagCSVOnline2", bTagCSVOnline2_, "bTagCSVOnline2/F");
//  outTree_->Branch("bTagCSVOnline3", bTagCSVOnline3_, "bTagCSVOnline3/F");
//  outTree_->Branch("bTagCSVOnline4", bTagCSVOnline4_, "bTagCSVOnline4/F");
//  outTree_->Branch("bTagCSVOffline1", bTagCSVOffline1_, "bTagCSVOffline1/F");
//  outTree_->Branch("bTagCSVOffline2", bTagCSVOffline2_, "bTagCSVOffline2/F");
//  outTree_->Branch("bTagCSVOffline3", bTagCSVOffline3_, "bTagCSVOffline3/F");
//  outTree_->Branch("bTagCSVOffline4", bTagCSVOffline4_, "bTagCSVOffline4/F");
  outTree_->Branch("nbtagCSVOnlinejets", &nbtagCSVOnlinejets_, "nbtagCSVOnlinejets/I");
  outTree_->Branch("bTagCSVOnline",  bTagCSVOnline_, "bTagCSVOnline[nbtagCSVOnlinejets]/F"); 
  outTree_->Branch("nbtagCSVOfflinejets", &nbtagCSVOfflinejets_, "nbtagCSVOfflinejets/I");
  outTree_->Branch("bTagCSVOffline",  bTagCSVOffline_, "bTagCSVOffline[nbtagCSVOfflinejets]/F");
//  outTree_->Branch("pfjetpt1", pfjetpt1_, "pfjetpt1/F");
//  outTree_->Branch("pfjetpt2", pfjetpt2_, "pfjetpt2/F");
//  outTree_->Branch("pfjetpt3", pfjetpt3_, "pfjetpt3/F");
//  outTree_->Branch("pfjetpt4", pfjetpt4_, "pfjetpt4/F");  
//  outTree_->Branch("pfjeteta1", pfjeteta1_, "pfjeteta1/F");
//  outTree_->Branch("pfjeteta2", pfjeteta2_, "pfjeteta2/F");
//  outTree_->Branch("pfjeteta3", pfjeteta3_, "pfjeteta3/F");
//  outTree_->Branch("pfjeteta4", pfjeteta4_, "pfjeteta4/F");
//  outTree_->Branch("pfjetphi1", pfjetphi1_, "pfjetphi1/F");
//  outTree_->Branch("pfjetphi2", pfjetphi2_, "pfjetphi2/F");
//  outTree_->Branch("pfjetphi3", pfjetphi3_, "pfjetphi3/F");
//  outTree_->Branch("pfjetphi4", pfjetphi4_, "pfjetphi4/F");
//  outTree_->Branch("pfjetmass1", pfjetmass1_, "pfjetmass1/F");
//  outTree_->Branch("pfjetmass2", pfjetmass2_, "pfjetmass2/F");
//  outTree_->Branch("pfjetmass3", pfjetmass3_, "pfjetmass3/F");
//  outTree_->Branch("pfjetmass4", pfjetmass4_, "pfjetmass4/F");
  outTree_->Branch("isptOrdered", &isptOrdered_, "isptOrdered/I");
  outTree_->Branch("npfjets",  &npfjets_,  "npfjets/I");
  outTree_->Branch("pfjetpt",  pfjetpt_,   "pfjetpt[npfjets]/F");
  outTree_->Branch("pfjeteta", pfjeteta_,   "pfjeteta[npfjets]/F");
  outTree_->Branch("pfjetphi", pfjetphi_,  "pfjetphi[npfjets]/F");
  outTree_->Branch("pfjetmass", pfjetmass_, "pfjetmass[npfjets]/F");
//  outTree_->Branch("nhltjets",   &nhltjets_,  "nhltjets/I");
//  outTree_->Branch("hltjetpt",   hltjetpt_,   "hltjetpt[nhltjets]/F");
//  outTree_->Branch("hltjeteta",  hltjeteta_,  "hltjeteta[nhltjets]/F");
//  outTree_->Branch("hltjetphi",  hltjetphi_,  "hltjetphi[nhltjets]/F");
//  outTree_->Branch("hltjetmass", hltjetmass_, "hltjetmass[nhltjets]/F");
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

//  edm::Handle<trigger::TriggerEvent> trigsummary;
//  iEvent.getByToken(trigsummaryToken_, trigsummary);

   // filtering by trigger path for the online objects
   // correct trigger path ensures objects were created => exist 
   // trigger::size_type trigFilterJet_Index       = trigsummary->filterIndex( trigFilterJet_ );
   // trigger::TriggerObjectCollection triggerObjects = trigsummary->getObjects();

//  if(trigFilterJet_Index < trigsummary->sizeFilters()) {
//    const trigger::Keys& trigKeys = trigsummary->filterKeys(trigFilterJet_Index); 
//    nhltjets_ = trigKeys.size();
//    // need kinematics for up to 4 leading jets
//    for(unsigned int ik = 0; ik < 4 && ik < trigKeys.size(); ++ik){ 
//      const trigger::TriggerObject& obj = triggerObjects[trigKeys[ik]];
//      hltjetpt_[ik] = obj.pt();
//      hltjeteta_[ik] = obj.eta();
//      hltjetphi_[ik] = obj.phi();
//      hltjetmass_[ik] = obj.mass();
//    }
//  }

  // kinematics from PFJets
  edm::Handle<reco::PFJetCollection> pfjetHandle;
  iEvent.getByToken(pfjetsToken_, pfjetHandle);
  const reco::PFJetCollection & pfjets = *(pfjetHandle.product());

  // want kinematics for up to 4 leading jets
  //if (pfjets.size() > 4) {
  //  npfjets_ = 4; 
  //}
  //else {
    npfjets_ = pfjets.size();
  //}
  
  // loop over jets
  for(int i = 0; i != npfjets_; ++i) {
    pfjetpt_[i]   = pfjets[i].pt();
    pfjeteta_[i]  = pfjets[i].eta();
    pfjetphi_[i]  = pfjets[i].phi();
    pfjetmass_[i] = pfjets[i].mass();
  }

  // check if pts are ordered
  // if (npfjets_ > 1) {
  //  for (int i = 1; i != npfjets_; ++i) {
  //    if (pfjetpt_[i] > pfjetpt_[i-1]) {
  //      isptOrdered_ = 0; 
  //    }
  //  }
  //}

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

  // want discriminants for up to 4 leading jets  
  //if (btagsCSVOffline.size() > 4) {
  //  nbtagCSVOfflinejets_ = 4; 
  //}
  //else {
    nbtagCSVOfflinejets_ = btagsCSVOffline.size();
  //}

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
