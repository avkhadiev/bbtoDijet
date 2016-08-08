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
//TGraphAsymmErrors* makeEffGraph(TH1F* pass, TH1F* total, bool debug=false) {
//
//  // make sure <pass> and <total> have the same binning!
//  int npoints = total->GetNbinsX();
//  float x[npoints], y[npoints], errx[npoints], erryl[npoints], erryh[npoints];
//  float npass = 0.0;
//  float ntotal = 0.0;
//  for(int ibin = 1; ibin < npoints+1; ibin++) {
//    x[ibin-1] = total->GetBinCenter(ibin);
//    npass = pass->GetBinContent(ibin);
//    ntotal = total->GetBinContent(ibin);
//    y[ibin-1] = ntotal < 1.0 ? 0.0 : npass/ntotal; // signal efficiency zero if no events in the bin
//    errx[ibin-1] = 0.0; // uncertainty in parameters is zero
//    if(y[ibin-1] == 0.0) {
//    } else {
//      if(debug) printf("npass = %3.1f, ntotal = %3.1f, eff = %4.2f", npass, ntotal, y[ibin-1]);
//      // the thid argument in ClopperPearson is the confidence interval; 0.683 is the value recommended by the PDG
//      // see https://root.cern.ch/doc/master/classTEfficiency.html#ae80c3189bac22b7ad15f57a1476ef75b
//      // fourth parameter is a bool: is the value for the upper bound?
//      erryl[ibin-1] = y[ibin-1] - TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, false); 
//      erryh[ibin-1] = TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, true) - y[ibin-1]; 
//    }
//  }
//
//  TGraphAsymmErrors *gr = new TGraphAsymmErrors(npoints, x, y, errx, errx, erryl, erryh);
//
//  return gr;
//
//}

void runEffPlots(std::vector<TString> infileNames)
{

  const char * tagName   = "HLT_DoubleJetsC100_p014_DoublePFJetsC100MaxDeta1p6_v2";
  const char * probeName = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v2"; 
  TChain *chain = new TChain("bbtoDijet/efficiencyTree");

  // add input file names
  for(const auto fname : infileNames) {
    chain->Add(fname);
  }

  int nevents = chain->GetEntries();
  printf("Will process %d events\n",nevents);

  // array sizes
  const size_t kMaxTriggerPass     = 1000    ;
  //const size_t kMaxCaloJet         = 1000    ;
  //const size_t kMaxPFJet           = 1000    ;
  const size_t kMaxBTagCSVOnline   = 1000    ;
  const size_t kMaxBTagCSVOffline  = 1000    ;

  // set of variables for 
  // trigger-related information
  unsigned int passTrigger                             = 0      ;    chain->SetBranchAddress(probeName,                    &passTrigger)                  ;
  unsigned int passControl                             = 0      ;    chain->SetBranchAddress(tagName,                      &passControl)                  ;
  // jet-related information 
  // CaloJet 
  //int          caloJetCounter                          = 0      ;    chain->SetBranchAddress("caloJetCounter",           &caloJetCounter)               ;
  //float        caloJetEta[kMaxCaloJet]                          ;    chain->SetBranchAddress("caloJetEta",               &caloJetEta)                   ; 
  //float        caloJetPhi[kMaxCaloJet]                          ;    chain->SetBranchAddress("caloJetPhi",               &caloJetPhi)                   ;  
  //double       caloJetEt[kMaxCaloJet]                           ;    chain->SetBranchAddress("caloJetEt",                &caloJetEt)                    ;  
  //// PFJet
  //int          pfJetCounter                            = 0      ;    chain->SetBranchAddress("pfJetCounter",             &pfJetCounter)                 ;
  //double       pfJetPt[kMaxPFJet]                               ;    chain->SetBranchAddress("pfJetPt",                  &pfJetPt)                      ;  
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
  // online 
  TEfficiency * effOn    = new TEfficiency("effOn",   "Efficiency(max, submax online CSV);Max online CSV;submax online CSV", 50, 0, 1.0, 50, 0, 1.0);
  TH2F * passHistoOn     = new TH2F("passHistoOn",    "Passing p84;max online CSV;submax oline CSV", 50, 0, 1.0, 50, 0, 1.0);
  TH2F * controlHistoOn  = new TH2F("controlHistoOn", "Control p84;max online CSV;submax online CSV", 50, 0, 1.0, 50, 0, 1.0); 
  // offline
  TEfficiency * effOff   = new TEfficiency("effOff",   "Efficiency(offline CSV matching with max, submax Online CSV);offline CSV matching with max online CSV;offline CSV matching with submax online CSV", 50, 0, 1.0, 50, 0, 1.0);
  TH2F * passHistoOff    = new TH2F("passHistoOff",    "Passing p84;offline CSV matching with max online CSV;offline CSV matching with submax online CSV", 50, 0, 1.0, 50, 0, 1.0);
  TH2F * controlHistoOff = new TH2F("controlHistoOff", "Control p84;offline CSV matching with max online CSV;offline CSV matching with submax online CSV", 50, 0, 1.0, 50, 0, 1.0);
  // efficiencies
  TH2D * corrHisto1      = new TH2D("corrHisto1",  "Offline vs. Online CSV [Max Online];max online CSV;matching offline CSV",       50, 0, 1, 50, 0, 1);
  TH2D * corrHisto2      = new TH2D("corrHisto2",  "Offline vs. Online CSV [Submax Online];submax online CSV;matching offline CSV", 50, 0, 1, 50, 0, 1);
  TH2D * corrHistoPt     = new TH2D("corrHistoPt", "Submax Online vs. Max Online CSV [p_{T} index];max online CSV p_{T} index;submax online CSV p_{T} index", 5, 0, 5, 5, 0, 5);
  
  // Loop over all jets in the event
  for(int ievent = 0; ievent < nevents; ievent++) {
    chain->GetEntry(ievent);
 
    // want to only consider up to 6 highest Et jets to match online selection
    int njets;
    if(bTagCSVOnlineJetCounter > 6){
    	njets = 6;
    }
    else{
        njets = bTagCSVOnlineJetCounter;
    }

    // need at least 2 b-tagged jets (other events passing the control have discriminant -10.0 and are not plotted)
    if((njets > 1) && passControl){
	
        int iOrderedEt[bTagCSVOnlineJetCounter];
        TMath::Sort(bTagCSVOnlineJetCounter, bTagCSVOnlineJetEt, iOrderedEt); // output array of indices corresponding to the decreasing order of et 
        int iOrderedPt[bTagCSVOnlineJetCounter];
        TMath::Sort(bTagCSVOnlineJetCounter, bTagCSVOnlineJetPt, iOrderedPt); // output array of indices corresponding to the decreasing order of pt 
        
        float minDeltaR[2]        = {10.0, 10.0}    ;
        float maxCSVOnline[2]     = {-10.0, -10.0}  ;
        int   iMaxCSVOnline[2]    = {0, 0}          ;
        float matchCSVOffline[2]  = {-10.0, -10.0}  ;
        int   iMatchCSVOffline[2] = {0, 0}          ;   
        
        // find max and submax online discriminants
        // CSV tagger considers jets with Et > 30, CSV matcher considers jets with Et > 80. 
        // to correct for bad input tag difference between CSV tagger and CSV matcher, enforce offline selection to 
        // emulate the orrect online behaviour 
        for(int ijet = 0; ijet < njets; ijet++) { 
            if(bTagCSVOnlineJetEt[iOrderedEt[ijet]] > 80.0){
                if(maxCSVOnline[0] < bTagCSVOnline[iOrderedEt[ijet]]){
                    iMaxCSVOnline[1] = iMaxCSVOnline[0]                ;
                    iMaxCSVOnline[0] = iOrderedEt[ijet]                ;
                    maxCSVOnline[1]  = maxCSVOnline[0]                 ;
                    maxCSVOnline[0]  = bTagCSVOnline[iOrderedEt[ijet]] ;
                }
                else if(maxCSVOnline[1] < bTagCSVOnline[iOrderedEt[ijet]]){
                    iMaxCSVOnline[1] = iOrderedEt[ijet]                            ;
                    maxCSVOnline[1]  = bTagCSVOnline[iOrderedEt[ijet]]             ;
                }
            }
            else { break; }
        }

	// match offline discriminants to max and submax online discriminants -- if both had Et > 80.0 GeV
        if(maxCSVOnline[1] > -10.0){ 
	    for(int itag = 0; itag < 2; itag ++) {
	        minDeltaR[itag] = 10.0 ;
	        for(int imatch = 0; imatch < bTagCSVOfflineJetCounter; imatch ++) {
	    	float iDeltaR = deltaR(bTagCSVOnlineJetEta[iMaxCSVOnline[itag]], bTagCSVOnlineJetPhi[iMaxCSVOnline[itag]], bTagCSVOfflineJetEta[imatch], bTagCSVOfflineJetPhi[imatch]) ;
	    	if(iDeltaR < minDeltaR[itag]){
	    	    minDeltaR[itag] = iDeltaR;
	    	    iMatchCSVOffline[itag] = imatch ; 
	    	    matchCSVOffline[itag]  = bTagCSVOffline[imatch] ;
	    	}
	        }
	    }
		
            // debugging 
            if(ievent % 1000 == 0) { 
               printf("********************************************* \n");
	       printf("Max online = %4.2f, Submax online = %4.2f \n", maxCSVOnline[0], maxCSVOnline[1]);
               printf("index in decreasing Et: iMax = %d, iSubmax = %d \n", iMaxCSVOnline[0], iMaxCSVOnline[1]);  
               printf("corresponding Et: Max = %4.2f, Submax = %4.2f \n", bTagCSVOnlineJetEt[iMaxCSVOnline[0]], bTagCSVOnlineJetEt[iMaxCSVOnline[1]]);
               printf("index in decreasing Pt: iMax = %d, iSubmax = %d \n", iOrderedPt[iMaxCSVOnline[0]], iOrderedPt[iMaxCSVOnline[1]]);  
               printf("corresponding Pt: Max = %4.2f, Submax = %4.2f \n", bTagCSVOnlineJetPt[iOrderedPt[iMaxCSVOnline[0]]], bTagCSVOnlineJetPt[iOrderedPt[iMaxCSVOnline[1]]]);
	       printf("Offline matching to max online = %4.2f, offline matching to submax online = %4.2f \n", matchCSVOffline[0], matchCSVOffline[1]);
	       printf("...with deltaR %f, %f, respectively \n", minDeltaR[0], minDeltaR[1]);
	    }

	    // fill the histograms if online jets were well-matched with offline jets
	    if((minDeltaR[0] < 0.4) && (minDeltaR[1] < 0.4)){
                if(ievent % 1000 == 0) {
	           printf("Filling histograms and plots \n");
                }
                // efficiency plots
	        effOn->Fill(passTrigger, maxCSVOnline[0], maxCSVOnline[1]);
	        effOff->Fill(passTrigger, matchCSVOffline[0], matchCSVOffline[1]);
 	        // histograms for events passing the control trigger
	        controlHistoOn->Fill(maxCSVOnline[0], maxCSVOnline[1]);	      
	        controlHistoOff->Fill(matchCSVOffline[0], matchCSVOffline[1]);
	        // offline-online CSV correlation histograms
	        corrHisto1->Fill(maxCSVOnline[0], matchCSVOffline[0]) ;	      
	        corrHisto2->Fill(maxCSVOnline[1], matchCSVOffline[1]) ;	      
	        // histograms for events passing the tagging trigger
                if(passTrigger){
	           passHistoOn->Fill(maxCSVOnline[0], maxCSVOnline[1]);
	           passHistoOff->Fill(matchCSVOffline[0], matchCSVOffline[1]);
                   // pT-index correlation histogram`
	           corrHistoPt->Fill(iOrderedPt[iMaxCSVOnline[0]], iOrderedPt[iMaxCSVOnline[1]]);
	        }
            }
	}
     }
  }

  // make plots
  // write histograms as tree branches
  printf("Writing histograms to the file \n"); 
  TFile outputFile("plots_CSV.root", "recreate");
  // efficiencies
  outputFile.mkdir("efficiencies");
  outputFile.cd("efficiencies");
  effOn->Write();
  effOff->Write();
  controlHistoOn->Write();
  controlHistoOff->Write();
  passHistoOn->Write();
  passHistoOff->Write();
  outputFile.cd();
  // correlations 
  outputFile.mkdir("correlations");
  outputFile.cd("correlations");
  corrHisto1->Write();
  corrHisto2->Write();
  corrHistoPt->Write();
  outputFile.cd(); 
  outputFile.Close();

  // writing out pdfs
  printf("Saving canvases as pdfs \n"); 
  gROOT->SetBatch()                               ; 
  gStyle->SetOptStat(0)                           ;
  TCanvas *c = MakeCanvas("c","",800,600)         ;
  // efficiencies
  controlHistoOn->GetXaxis()->SetRangeUser(0.5, 1)  ;
  controlHistoOn->GetYaxis()->SetRangeUser(0.5, 1)  ;
  controlHistoOn->Draw("COLZ")                    ;
  c->SaveAs("histo_on_control_p84.pdf")           ;
  controlHistoOff->GetXaxis()->SetRangeUser(0.5, 1) ;
  controlHistoOff->GetYaxis()->SetRangeUser(0.5, 1) ;
  controlHistoOff->Draw("COLZ")                   ;
  c->SaveAs("histo_off_control_p84.pdf")          ;
  passHistoOn->GetXaxis()->SetRangeUser(0.83, 1)  ;
  passHistoOn->GetYaxis()->SetRangeUser(0.83, 1)  ;
  passHistoOn->Draw("COLZ")                       ;
  c->SaveAs("histo_on_tagging_p84.pdf")           ;
  passHistoOff->GetXaxis()->SetRangeUser(0.83, 1) ;
  passHistoOff->GetYaxis()->SetRangeUser(0.83, 1) ;
  passHistoOff->Draw("COLZ")                      ;
  c->SaveAs("histo_off_tagging_p84.pdf")          ;
  effOn->Draw("COLZ")                             ;
  c->SaveAs("eff_on_CSVp84.pdf")                  ;
  effOff->Draw("COLZ")                             ;
  c->SaveAs("eff_off_CSVp84.pdf")                  ;
  // correlations
  int zMin                                        ;
  zMin = ceil(0.0005 * corrHisto1->GetEntries())  ; 
  corrHisto1->SetMinimum(zMin)                    ;
  corrHisto1->Draw("COLZ")                        ;
  c->SaveAs("corr_1_CSVp84.pdf")                  ;
  zMin = ceil(0.00005  * corrHisto2->GetEntries()); 
  corrHisto2->SetMinimum(zMin)                    ;
  corrHisto2->Draw("COLZ")                        ;
  c->SaveAs("corr_2_CSVp84.pdf")                  ;
  c->SetGrid(5, 5)                                ;
  gStyle->SetOptStat()                            ;
  corrHistoPt->SetMarkerSize(1.2)                 ;
  corrHistoPt->SetMarkerColor(2)                  ;
  corrHistoPt->Draw("TEXT")                       ;
  c->SaveAs("corr_ptIndex_CSVp84.pdf")            ;
}

void makeEffPlots()
{

  std::vector<TString> filelist;

  // "root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root"
  // specify output tree
  filelist.push_back("root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root");

  runEffPlots(filelist);
}
