#include <string>

#include "TFile.h"
#include "TMath.h"

void sortTree() {
    TFile f1("test_bTagDijetV11.root");
    TTree *tree = (TTree*)f1.Get("bbtoDijet/efficiencyTree");
    Int_t nentries = (Int_t)tree->GetEntries();
    // drawing variable param with no graphics option.
    // variable param stored in array fV1 (see TTree::Draw)
    tree->Draw("caloJetEt","","goff");
    Int_t *index = new Int_t[nentries];
    // sort array containing param in decreasing order
    // the array index contains the entry numbers in decreasing order of param 
    TMath::Sort(nentries, tree->GetV1(), index);
    // open new file to store the sorted Tree
    TFile f2("sorted.root","recreate");
    // create an empty clone of the original tree
    TTree *tsorted = (TTree*)tree->CloneTree(0);
    for (Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(index[i]);
        tsorted->Fill();
    }
    tsorted->Write();
    delete [] index;
}                 
