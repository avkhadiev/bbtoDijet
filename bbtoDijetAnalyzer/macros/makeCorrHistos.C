#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <vector>
#include <string>
#include <math.h>       /* ceil */
#include "Style.hh"
#endif

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
 float deltaPhi = TMath::Abs(phi1-phi2);
 float deltaEta = eta1-eta2;
 if(deltaPhi > TMath::Pi())
   deltaPhi = TMath::TwoPi() - deltaPhi;
 return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

// Makes a turn-on curve once the histograms (with same binnings!) are ready
// Error calculation: http://indico.cern.ch/event/109265/contributions/1312303/attachments/29409/42521/trigeff-20101019.pdf
// also see https://root.cern.ch/doc/master/classTEfficiency.html
TGraphAsymmErrors* makeEffGraph(TH1F* pass, TH1F* total, bool debug=false) {

  // make sure <pass> and <total> have the same binning!
  int npoints = total->GetNbinsX();
  float x[npoints], y[npoints], errx[npoints], erryl[npoints], erryh[npoints];
  float npass = 0.0;
  float ntotal = 0.0;
  for(int ibin = 1; ibin < npoints+1; ibin++) {
    x[ibin-1] = total->GetBinCenter(ibin);
    npass = pass->GetBinContent(ibin);
    ntotal = total->GetBinContent(ibin);
    y[ibin-1] = ntotal < 1.0 ? 0.0 : npass/ntotal; // signal efficiency zero if no events in the bin
    errx[ibin-1] = 0.0; // uncertainty in parameters is zero
    if(y[ibin-1] == 0.0) {
    } else {
      if(debug) printf("npass = %3.1f, ntotal = %3.1f, eff = %4.2f", npass, ntotal, y[ibin-1]);
      // the thid argument in ClopperPearson is the confidence interval; 0.683 is the value recommended by the PDG
      // see https://root.cern.ch/doc/master/classTEfficiency.html#ae80c3189bac22b7ad15f57a1476ef75b
      // fourth parameter is a bool: is the value for the upper bound?
      erryl[ibin-1] = y[ibin-1] - TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, false); 
      erryh[ibin-1] = TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, true) - y[ibin-1]; 
    }
  }

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(npoints, x, y, errx, errx, erryl, erryh);

  return gr;

}

void runCorrHistos(std::vector<TString> infileNames)
{

  const char * tagName   = "HLT_DoubleJetsC100_p014_DoublePFJetsC100MaxDeta1p6_v2" ;// "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v2";
  const char * probeName = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v2" ;//"HLT_DoubleJetsC100_p026_DoublePFJetsC160_v2"; 
  TChain *chain = new TChain("bbtoDijet/efficiencyTree");

  // add input file names
  for(const auto fname : infileNames) {
    chain->Add(fname);
  }

  int nevents = chain->GetEntries();
  printf("Will process %d events\n",nevents);

  // array sizes
  const size_t kMaxTriggerPass     = 1000    ;
  const size_t kMaxCaloJet         = 1000    ;
  const size_t kMaxPFJet           = 1000    ;
  const size_t kMaxBTagCSVOnline   = 1000    ;
  const size_t kMaxBTagCSVOffline  = 1000    ;

  // set of variables for 
  // trigger-related information
  unsigned int passTrigger                             = 0      ;    chain->SetBranchAddress(probeName,                    &passTrigger)                ;
  unsigned int passControl                             = 0      ;    chain->SetBranchAddress(tagName,                      &passControl)                ;
  // jet-related information 
  // CaloJet 
  int          caloJetCounter                          = 0      ;    chain->SetBranchAddress("caloJetCounter",           &caloJetCounter)               ;
  float        caloJetEta[kMaxCaloJet]                          ;    chain->SetBranchAddress("caloJetEta",               &caloJetEta)                   ; 
  float        caloJetPhi[kMaxCaloJet]                          ;    chain->SetBranchAddress("caloJetPhi",               &caloJetPhi)                   ;  
  double       caloJetEt[kMaxCaloJet]                           ;    chain->SetBranchAddress("caloJetEt",                &caloJetEt)                    ;  
  // PFJet
  int          pfJetCounter                            = 0      ;    chain->SetBranchAddress("pfJetCounter",             &pfJetCounter)                 ;
  double       pfJetPt[kMaxPFJet]                               ;    chain->SetBranchAddress("pfJetPt",                  &pfJetPt)                      ;  
  // b-jet-tagging information 
  // CSVOnline
  int          bTagCSVOnlineJetCounter                 = 0      ;    chain->SetBranchAddress("bTagCSVOnlineJetCounter",  &bTagCSVOnlineJetCounter)      ;
  float        bTagCSVOnlineJetEta[kMaxBTagCSVOnline]           ;    chain->SetBranchAddress("bTagCSVOnlineJetEta",      &bTagCSVOnlineJetEta)          ; 
  float        bTagCSVOnlineJetPhi[kMaxBTagCSVOnline]           ;    chain->SetBranchAddress("bTagCSVOnlineJetPhi",      &bTagCSVOnlineJetPhi)          ; 
  double       bTagCSVOnlineJetEt[kMaxBTagCSVOnline]            ;    chain->SetBranchAddress("bTagCSVOnlineJetEt",       &bTagCSVOnlineJetEt)           ;  
  double       bTagCSVOnlineJetPt[kMaxBTagCSVOnline]            ;    chain->SetBranchAddress("bTagCSVOnlineJetPt",       &bTagCSVOnlineJetPt)           ;  
  float        bTagCSVOnline[kMaxBTagCSVOnline]                 ;    chain->SetBranchAddress("bTagCSVOnline",            &bTagCSVOnline)                ;
  // CSVOffline
  int          bTagCSVOfflineJetCounter                = 0      ;    chain->SetBranchAddress("bTagCSVOfflineJetCounter", &bTagCSVOfflineJetCounter)     ;
  float        bTagCSVOfflineJetEta[kMaxBTagCSVOffline]         ;    chain->SetBranchAddress("bTagCSVOfflineJetEta",     &bTagCSVOfflineJetEta)         ; 
  float        bTagCSVOfflineJetPhi[kMaxBTagCSVOffline]         ;    chain->SetBranchAddress("bTagCSVOfflineJetPhi",     &bTagCSVOfflineJetPhi)         ; 
  double       bTagCSVOfflineJetEt[kMaxBTagCSVOffline]          ;    chain->SetBranchAddress("bTagCSVOfflineJetEt",      &bTagCSVOfflineJetEt)          ;  
  double       bTagCSVOfflineJetPt[kMaxBTagCSVOffline]          ;    chain->SetBranchAddress("bTagCSVOfflineJetPt",      &bTagCSVOfflineJetPt)          ;  
  float        bTagCSVOffline[kMaxBTagCSVOffline]               ;    chain->SetBranchAddress("bTagCSVOffline",           &bTagCSVOffline)               ;

  // declare histograms
  TH2D * correlationHisto1  = new TH2D("correlationHisto1",  "Offline vs. Online CSV [Max];Online CSV;Offline CSV;Accepted",    50, 0, 1, 50, 0, 1);
  TH2D * correlationHisto2  = new TH2D("correlationHisto2",  "Offline vs. Online CSV [SubMax];Online CSV;Offline CSV;Accepted", 50, 0, 1, 50, 0, 1);
  TH2D * correlationHistoPt = new TH2D("correlationHistoPt", "Submax Online vs. Max Online CSV [p_{T} index];Max Online CSV p_{T} index;Submax Online CSV p_{T} index;Accepted", 6, 0, 5, 6, 0, 5);
  // Loop over all jets in the event
  for(int ievent = 0; ievent < nevents; ievent++) {
    chain->GetEntry(ievent);
 
  //  loop over all jets in the event, store indices of jets with up to 8 highest Et values
  //  double maxEt[8]  = {0.0}
  //  int    iMaxEt[8] = {0}
  //  for(int ijet = 0; ijet < caloJetCounter; ijet++) {
  //      for(int iEt = 0; iEt < 8; iEt++) {
  //          if(maxEt[iEt] < caloJetEt[ijet]) {
  //              for(iShift = 7; iShift > iEt; iShift--) {
  //                  iMaxEt[iShift] = iMaxEt[iShift - 1]
  //                  maxEt[iShift]  = maxEt[iShift - 1]
  //              }
  //              iMaxEt[iEt] = ijet
  //              maxEt[iEt]  = caloJetEt[ijet]
  //              break;
  //          }
  //      }
        
  // ****************************************************************************************
    // filter bTagCSVOnline in decreasing order of bTagCSVOnlineJetPt -- it's not ordered originally!
    int iOrderedPt[bTagCSVOnlineJetCounter];
    TMath::Sort(bTagCSVOnlineJetCounter, bTagCSVOnlineJetPt, iOrderedPt); // output array of indices corresponding to the decreasing order of pt 

    if(passControl) {
	float minDeltaR[2]        = {10.0, 10.0}    ;
        float maxCSVOnline[2]     = {-10.0, -10.0}  ;
        int   iMaxCSVOnline[2]    = {0, 0}          ;
        float matchCSVOffline[2]  = {-10.0, -10.0}  ;
        int   iMatchCSVOffline[2] = {0, 0}          ;
        // find max and submax online discriminants
        for(int ijet = 0; ijet < bTagCSVOnlineJetCounter; ijet++) { 
            if(maxCSVOnline[0] < bTagCSVOnline[iOrderedPt[ijet]]) {
                iMaxCSVOnline[1] = iMaxCSVOnline[0]    ;
                iMaxCSVOnline[0] = iOrderedPt[ijet]                ;
                maxCSVOnline[1]  = maxCSVOnline[0]     ;
                maxCSVOnline[0]  = bTagCSVOnline[iOrderedPt[ijet]] ;
            }
            else if(maxCSVOnline[1] < bTagCSVOnline[iOrderedPt[ijet]]){
                iMaxCSVOnline[1] = iOrderedPt[ijet]                ;
                maxCSVOnline[1]  = bTagCSVOnline[iOrderedPt[ijet]] ;
            }
        }
        // match offline discriminants to max and submax online discriminants
        for(int itag = 0; itag < 2; itag ++) {
            minDeltaR[itag] = 10.0 ;
            for(int imatch = 0; imatch < bTagCSVOfflineJetCounter; imatch ++) {
                float iDeltaR = deltaR(bTagCSVOnlineJetEta[iMaxCSVOnline[itag]], bTagCSVOnlineJetPhi[iMaxCSVOnline[itag]], bTagCSVOfflineJetEta[imatch], bTagCSVOfflineJetPhi[imatch]) ;
                if(iDeltaR < minDeltaR[itag]){
                    minDeltaR[itag] = iDeltaR;
                    matchCSVOffline[itag]  = bTagCSVOffline[imatch] ;
                    iMatchCSVOffline[itag] = imatch ; 
                }
            }
	        // debugging 
            if(ievent % 100 == 0) { 
	    	printf("minDeltaR = %f \n", minDeltaR[itag]);
            printf("Online = %4.2f, Offline = %4.2f \n", maxCSVOnline[itag], matchCSVOffline[itag]);
        	printf("Max online pt index = %d, Submax online pt index = %d \n", iMaxCSVOnline[0], iMaxCSVOnline[1]);
        	printf("Max online = %4.2f, Submax online = %4.2f \n", maxCSVOnline[0], maxCSVOnline[1]);
        	printf("Pt for max online = %4.2f, Pt for submax online = %4.2f \n", bTagCSVOnlineJetPt[iMaxCSVOnline[0]],bTagCSVOnlineJetPt[iMaxCSVOnline[1]]);
            }
	}
	if(minDeltaR[0] < 0.1) { 
        correlationHisto1->Fill(maxCSVOnline[0], matchCSVOffline[0]) ;
        if(passTrigger){
        correlationHistoPt->Fill(iMaxCSVOnline[0], iMaxCSVOnline[1]) ;
        }
    }
    if(minDeltaR[1] < 0.1) { correlationHisto2->Fill(maxCSVOnline[1], matchCSVOffline[1]) ; }
    }
  }

  // make plots
  // write histograms as tree branches 
  TFile outputFile("corrHistos.root", "recreate")                       ;
  outputFile.mkdir("histos_");
  outputFile.cd("histos_");
  correlationHisto1->Write();
  correlationHisto2->Write();
  correlationHistoPt->Write();
  outputFile.cd();
  outputFile.Close();

  //gStyle->SetOptStat(0)                           ;
  //int zMin, zMax                                  ; 
  //TCanvas *c = MakeCanvas("corrCanvas","",800,600);
  //zMin = ceil(0.0001 * correlationHisto1->GetEntries()) ; 
  ////  zMax = ceil(0.1    * correlationHisto1->GetEntries()) ; 
  //correlationHisto1->SetMinimum(zMin)             ;
  //correlationHisto1->Draw("COLZ")                 ;
  //c->SaveAs("correlation_1.pdf")                  ;
  //zMin = ceil(0.0001 * correlationHisto2->GetEntries()) ; 
  //// zMax = ceil(0.1    * correlationHisto2->GetEntries()) ; 
  //correlationHisto2->Draw("COLZ")                 ;
  //correlationHisto2->SetMinimum(zMin)             ;
  //c->SaveAs("correlation_2.pdf")                  ;
}

void makeCorrHistos()
{

  std::vector<TString> filelist;

  // "root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root"
  // specify output tree
  filelist.push_back("root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root");

  runCorrHistos(filelist);

}
