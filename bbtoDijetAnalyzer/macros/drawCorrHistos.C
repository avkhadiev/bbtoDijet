#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2D.h"
#endif

void drawCorrHistos() {
    TFile *f = new TFile("corrHistos.root");
    TTree *T = (TTree*)f->Get("histosTree");
    TH2D *correlationHisto1 = 0;
    TH2D *correlationHisto2 = 0;
    T->SetBranchAddress("correlationHisto1", &correlationHisto1);
    T->SetBranchAddress("correlationHisto2", &correlationHisto2);
    TCanvas *c1 = new TCanvas("c1", "Offline vs. Online CSV [Max]", 800, 600);
    c1->SetFillColor( 19 );
    correlationHisto1->Draw("COLZ");
    c1->SaveAs("c1.pdf");
    TCanvas *c2 = new TCanvas("c2", "Offline vs. Online CSV [Submax]", 800, 600);
    c2->SetFillColor( 19 );
    correlationHisto2->Draw("COLZ");
    c2->SaveAs("c2.pdf");
}
