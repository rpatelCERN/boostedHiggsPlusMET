//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 11 12:42:38 2018 by ROOT version 6.06/09
// from TTree newtree/newtree
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef BoostedLoop_h
#define BoostedLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class BoostedLoop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           boostBool1;
   Int_t           boostBool2;
   Int_t           boostBool3;
   Int_t           boostBool4;
   Int_t           boostBool5;
   Int_t           boostBool6;
   Int_t           boostBool7;
   Int_t           nGenHs;
   Int_t           BTags;
   Float_t         HT;
   Float_t         MET;
   Double_t        Weight;


   // List of branches
   TBranch        *b_boostBool1;   //!
   TBranch        *b_boostBool2;   //!
   TBranch        *b_boostBool3;   //!
   TBranch        *b_boostBool4;   //!
   TBranch        *b_boostBool5;   //!
   TBranch        *b_boostBool6;   //!
   TBranch        *b_boostBool7;   //!
   TBranch        *b_nGenHs;   //!
   TBranch        *b_BTags;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_Weight;   //!


   BoostedLoop(TTree *tree=0);
   virtual ~BoostedLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual vector<float>    Loop(float EventWeight, string file,bool save_files);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BoostedLoop_cxx
BoostedLoop::BoostedLoop(TTree *tree) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output.root");
      }
      f->GetObject("QCD",tree);

   }
   Init(tree);
}

BoostedLoop::~BoostedLoop() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BoostedLoop::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BoostedLoop::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BoostedLoop::Init(TTree *tree) {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).


   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("boostBool1", &boostBool1, &b_boostBool1);
   fChain->SetBranchAddress("boostBool2", &boostBool2, &b_boostBool2);
   fChain->SetBranchAddress("boostBool3", &boostBool3, &b_boostBool3);
   fChain->SetBranchAddress("boostBool4", &boostBool4, &b_boostBool4);
   fChain->SetBranchAddress("boostBool5", &boostBool5, &b_boostBool5);
   fChain->SetBranchAddress("boostBool6", &boostBool6, &b_boostBool6);
   fChain->SetBranchAddress("boostBool7", &boostBool7, &b_boostBool7);
   fChain->SetBranchAddress("nGenHs", &nGenHs, &b_nGenHs);
   fChain->SetBranchAddress("BTags", &BTags, &b_BTags);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);


   Notify();
}

Bool_t BoostedLoop::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BoostedLoop::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BoostedLoop::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BoostedLoop_cxx
