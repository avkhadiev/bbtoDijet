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

  // declare histograms
  TH2D * distributionCSV = new TH2D("distributionCSV", "Control p84;max online CSV;submax online CSV", 50, 0, 1, 50, 0, 1);
      xBins = distribitionCSV->GetXaxis->GetNbins();
      yBins = distributionCSV->GetYaxis->GetNbins();
      printf("xBins = %d, yBins = %d \n", xBins, yBins);
  graphSize = xBins * yBins;
  printf("declaring graph with size %d \n", graphSize); 
  TGraph2D rateGraph  = new TGraph2D( graphSize);  

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
        
        float maxCSVOnline[2]     = {-10.0, -10.0}  ;
        int   iMaxCSVOnline[2]    = {0, 0}          ;
        
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
	
            // debugging 
            if(ievent % 1000 == 0) { 
               printf("********************************************* \n");
	       printf("Max online = %4.2f, Submax online = %4.2f \n", maxCSVOnline[0], maxCSVOnline[1]);
               printf("index in decreasing Et: iMax = %d, iSubmax = %d \n", iMaxCSVOnline[0], iMaxCSVOnline[1]);  
               printf("corresponding Et: Max = %4.2f, Submax = %4.2f \n", bTagCSVOnlineJetEt[iMaxCSVOnline[0]], bTagCSVOnlineJetEt[iMaxCSVOnline[1]]);
	    }

	    // fill the histograms if online jets were well-matched with offline jets
                if(ievent % 1000 == 0) {
	           printf("Filling histograms and plots \n");
		}
                distributionCSV->Fill(MaxCSVOnline[0], MaxCSVOnline[1]);		
	}
     }
  }

  printf("CSV distribution built! \n Estimating rate now... \n");

  float xCoordinate, yCoordinate;
  int value;
  int i = 0;
  for(x = 0; x < xBins; x++){
      for(y = 0; y < yBins; y ++{
          xCoordinate = distribitionCSV->GetXaxis()->GetBinCenter(x);
          yCoordinate = distributionCSV->GetYaxis()->GetBinCenter(y);
          value = distibutionCSV->Integral(x, xBins-1, y, yBins - 1);
          // debugging 
          printf("xCoordinate is %4.2f, yCoordinate is %4.2f \n", xCoordinate, yCoordinate) ;
          printf("Calculated integral yields %d \n", value) ;
          rateGraph->SetPoint(i, xCoordinate, yCoordinate, value) ;
          i = ++i;
      }
  }

  // make plots
  // write histograms as tree branches 
  TFile outputFile("rate.root", "recreate");
  outputFile.mkdir("rates");
  outputFile.cd("rates");
  distributionCSV->Write();
  rateGraph->Write();
  outputFile.cd();
  outputFile.Close();
}

void estimateRate()
{

  std::vector<TString> filelist;

  // "root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root"
  // specify output tree
  filelist.push_back("root://cmsxrootd.fnal.gov//store/user/aavkhadi/JetHT/bbtoDijetV11/hlt_bTagDijetV11.root");

  makeGraph(filelist);

}
