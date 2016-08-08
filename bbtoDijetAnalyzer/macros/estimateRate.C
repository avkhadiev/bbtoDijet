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


void makeGraph(std::vector<TString> infileNames)
{

  const char * controlName   = "HLT_DoubleJetsC100_p014_DoublePFJetsC100MaxDeta1p6_v2";
  TChain *chain = new TChain("bbtoDijet/efficiencyTree");

  // add input file names
  for(const auto fname : infileNames) {
    chain->Add(fname);
  }

  int nevents = chain->GetEntries();
  printf("Will process %d events\n",nevents);

  // array sizes
  const size_t kMaxTriggerPass     = 1000    ;
  const size_t kMaxBTagCSVOnline   = 1000    ;

  // trigger-related information
  unsigned int passControl                             = 0      ;    chain->SetBranchAddress(controlName,                &passControl)                  ;
  // b-jet-tagging information 
  int          bTagCSVOnlineJetCounter                 = 0      ;    chain->SetBranchAddress("bTagCSVOnlineJetCounter",  &bTagCSVOnlineJetCounter)      ;
  double       bTagCSVOnlineJetEt[kMaxBTagCSVOnline]            ;    chain->SetBranchAddress("bTagCSVOnlineJetEt",       &bTagCSVOnlineJetEt)           ;  
  float        bTagCSVOnline[kMaxBTagCSVOnline]                 ;    chain->SetBranchAddress("bTagCSVOnline",            &bTagCSVOnline)                ;

  // declare histograms and graphs
  
  // Loop over all jets in the event
  for(int ievent = 0; ievent < nevents; ievent++) {
    chain->GetEntry(ievent);
 
    int iOrderedEt[bTagCSVOnlineJetCounter];
    TMath::Sort(bTagCSVOnlineJetCounter, bTagCSVOnlineJetEt, iOrderedEt); // output array of indices corresponding to the decreasing order of pt 
    int   njets                                 ;
    float maxCSVOnline[2]     = {-10.0, -10.0}  ;
    int   iMaxCSVOnline[2]    = {0, 0}          ;
    // want to only consider up to 6 highest Et jets to match online selection
    if(bTagCSVOnlineJetCounter > 6){
    	njets = 6;
    }
    else{
        njets = bTagCSVOnlineJetCounter;
    }
    // find max and submax online discriminants
    // CSV tagger considers jets with Et > 30, CSV matcher considers jets with Et > 80. 
    // to correct for bad input tag difference between CSV tagger and CSV matcher, enforce offline selection to 
    // emulate the orrect online behaviour 
    for(int ijet = 0; ijet < njets; ijet++) { 
        if(bTagCSVOnlineJetEt[iOrderedEt[ijet]] > 80.0){
            if(maxCSVOnline[0] < bTagCSVOnline[iOrderedEt[ijet]]) {
                iMaxCSVOnline[1] = iMaxCSVOnline[0]    ;
                iMaxCSVOnline[0] = iOrderedEt[ijet]                ;
                maxCSVOnline[1]  = maxCSVOnline[0]     ;
                maxCSVOnline[0]  = bTagCSVOnline[iOrderedEt[ijet]] ;
            }
            else if(maxCSVOnline[1] < bTagCSVOnline[iOrderedEt[ijet]]){
                iMaxCSVOnline[1] = iOrderedEt[ijet]                            ;
                maxCSVOnline[1]  = bTagCSVOnline[iOrderedEt[ijet]]             ;
            }
        }
	else { break; }
    }
    // debugging 
    if(ievent % 100 == 0) { 
    printf("Max online = %4.2f, Submax online = %4.2f \n", maxCSVOnline[0], maxCSVOnline[1]);
    printf("index in decreasing Et: iMax = %d, iSubmax = %d \n", iMaxCSVOnline[0], iMaxCSVOnline[1]);  
    printf("corresponding Et: Max = %4.2f, Submax = %4.2f \n", bTagCSVOnlineJetEt[iMaxCSVOnline[0]], bTagCSVOnlineJetEt[iMaxCSVOnline[1]]);
    }
    // current control (probe with CSV matcher disabled) allows jets with Et in (30.0, 80.0); 
    // to obtain adequate behavior, the search for max CSV was corrected to match the CSV filter selection
    if(passControl){
        effPlot->Fill(passTrigger, maxCSVOnline[0], maxCSVOnline[1]);
	controlHisto->Fill(maxCSVOnline[0], maxCSVOnline[1]);	      
	if(passTrigger){
	   passHisto->Fill(maxCSVOnline[0], maxCSVOnline[1]);
	}
	else {
	   missHisto->Fill(maxCSVOnline[0], maxCSVOnline[1]);
	}
    }
  }

  // make plots
  // write histograms as tree branches 
  TFile outputFile("plotEfficiencies.root", "recreate");
  outputFile.mkdir("efficiencies_");
  outputFile.cd("efficiencies_");
  effPlot->Write();
  controlHisto->Write();
  passHisto->Write();
  missHisto->Write();	  
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

void estimateRate()
{

  std::vector<TString> filelist;

  // "root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root"
  // specify output tree
  filelist.push_back("root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root");

  makeGraph(filelist);

}
