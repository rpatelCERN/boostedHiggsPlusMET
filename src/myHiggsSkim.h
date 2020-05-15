//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 23 07:37:49 2019 by ROOT version 6.06/09
// from TTree newtree/newtree
// found on file: output.root
//////////////////////////////////////////////////////////

#ifndef myHiggsSkim_h
#define myHiggsSkim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class myHiggsSkim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nAK4;
   Int_t           nDiAK8;
   Int_t           nBoostedH;
   Int_t           nResolvedH;
   Int_t           nGenZs;
   Int_t           nGenHs;
   Int_t           BTagsL;
   Int_t           BTagsM;
   Int_t           BTagsT;
   Float_t         H1_pt_res;
   Float_t         H2_pt_res;
   Float_t         H1_pt_boost;
   Float_t         H2_pt_boost;
   Float_t         MET;
   Float_t         HT;
   Float_t         res_J1btag;
   Float_t         res_J2btag;
   Float_t         res_J3btag;
   Float_t         res_J4btag;
   Float_t         res_avgM;
   Float_t         res_deltaM;
   Float_t         res_deltaRmax;
   Float_t         Weight;

   // List of branches
   TBranch        *b_nAK4;   //!
   TBranch        *b_nDiAK8;   //!
   TBranch        *b_nBoostedH;   //!
   TBranch        *b_nResolvedH;   //!
   TBranch        *b_nGenZs;   //!
   TBranch        *b_nGenHs;   //!
   TBranch        *b_BTagsL;   //!
   TBranch        *b_BTagsM;   //!
   TBranch        *b_BTagsT;   //!
   TBranch        *b_H1_pt_res;   //!
   TBranch        *b_H2_pt_res;   //!
   TBranch        *b_H1_pt_boost;   //!
   TBranch        *b_H2_pt_boost;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_res_J1btag;   //!
   TBranch        *b_res_J2btag;   //!
   TBranch        *b_res_J3btag;   //!
   TBranch        *b_res_J4btag;   //!
   TBranch        *b_res_avgM;   //!
   TBranch        *b_res_deltaM;   //!
   TBranch        *b_res_deltaRmax;   //!
   TBranch        *b_Weight;   //!

   myHiggsSkim(TTree *tree=0);
   virtual ~myHiggsSkim();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual vector<TH1F*> Loop(float EventWeight, TString file,bool save_files);
   virtual void eff(float EventWeight, TString process);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef myHiggsSkim_cxx
myHiggsSkim::myHiggsSkim(TTree *tree) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TString versionFile = "V18";
     TString thisFile = "output_"+version+".root";
     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(thisFile);
      if (!f || !f->IsOpen()) {
        f = new TFile(thisFile);
      }
      f->GetObject("newtree",tree);
   }
   Init(tree);
}

myHiggsSkim::~myHiggsSkim() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myHiggsSkim::GetEntry(Long64_t entry) { // Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t myHiggsSkim::LoadTree(Long64_t entry) { // Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void myHiggsSkim::Init(TTree *tree) {
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

   fChain->SetBranchAddress("nAK4", &nAK4, &b_nAK4);
   fChain->SetBranchAddress("nDiAK8", &nDiAK8, &b_nDiAK8);
   fChain->SetBranchAddress("nBoostedH", &nBoostedH, &b_nBoostedH);
   fChain->SetBranchAddress("nResolvedH", &nResolvedH, &b_nResolvedH);
   fChain->SetBranchAddress("nGenZs", &nGenZs, &b_nGenZs);
   fChain->SetBranchAddress("nGenHs", &nGenHs, &b_nGenHs);
   fChain->SetBranchAddress("BTagsL", &BTagsL, &b_BTagsL);
   fChain->SetBranchAddress("BTagsM", &BTagsM, &b_BTagsM);
   fChain->SetBranchAddress("BTagsT", &BTagsT, &b_BTagsT);
   fChain->SetBranchAddress("H1_pt_res", &H1_pt_res, &b_H1_pt_res);
   fChain->SetBranchAddress("H2_pt_res", &H2_pt_res, &b_H2_pt_res);
   fChain->SetBranchAddress("H1_pt_boost", &H1_pt_boost, &b_H1_pt_boost);
   fChain->SetBranchAddress("H2_pt_boost", &H2_pt_boost, &b_H2_pt_boost);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("res_J1btag", &res_J1btag, &b_res_J1btag);
   fChain->SetBranchAddress("res_J2btag", &res_J2btag, &b_res_J2btag);
   fChain->SetBranchAddress("res_J3btag", &res_J3btag, &b_res_J3btag);
   fChain->SetBranchAddress("res_J4btag", &res_J4btag, &b_res_J4btag);
   fChain->SetBranchAddress("res_avgM", &res_avgM, &b_res_avgM);
   fChain->SetBranchAddress("res_deltaM", &res_deltaM, &b_res_deltaM);
   fChain->SetBranchAddress("res_deltaRmax", &res_deltaRmax, &b_res_deltaRmax);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);

   Notify();
}

Bool_t myHiggsSkim::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myHiggsSkim::Show(Long64_t entry) { // Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t myHiggsSkim::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef myHiggsSkim_cxx
