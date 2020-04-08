//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 24 11:03:21 2020 by ROOT version 6.16/00
// from TTree BornTree/Born Tree
// found on file: BornTreeBrem.root
//////////////////////////////////////////////////////////

#ifndef BornAnalyzer_h
#define BornAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BornAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        x;
   Double_t        theta1;
   Double_t        theta2;
   Double_t        phi1;
   Double_t        phi2;
   Double_t        E0;
   Double_t        CrossSectionRatio;
   Double_t        HyperGeometricCrossSection;
   Int_t           Attempts;
   Int_t           ev;

   // List of branches
   TBranch        *b_x;   //!
   TBranch        *b_theta1;   //!
   TBranch        *b_theta2;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_E0;   //!
   TBranch        *b_CrossSectionRatio;   //!
   TBranch        *b_HyperGeometricCrossSection_t;   //!
   TBranch        *b_Attempts;   //!
   TBranch        *b_ev;   //!

   BornAnalyzer(TTree *tree=0);
   virtual ~BornAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BornAnalyzer_cxx
BornAnalyzer::BornAnalyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BornTreeBrem.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("BornTreeBrem.root");
      }
      f->GetObject("BornTree",tree);

   }
   Init(tree);
}

BornAnalyzer::~BornAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BornAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BornAnalyzer::LoadTree(Long64_t entry)
{
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

void BornAnalyzer::Init(TTree *tree)
{
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

   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("theta1", &theta1, &b_theta1);
   fChain->SetBranchAddress("theta2", &theta2, &b_theta2);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("E0", &E0, &b_E0);
   fChain->SetBranchAddress("CrossSectionRatio", &CrossSectionRatio, &b_CrossSectionRatio);
   fChain->SetBranchAddress("HyperGeometricCrossSection", &HyperGeometricCrossSection, &b_HyperGeometricCrossSection_t);
   fChain->SetBranchAddress("Attempts", &Attempts, &b_Attempts);
   fChain->SetBranchAddress("ev", &ev, &b_ev);
   Notify();
}

Bool_t BornAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BornAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BornAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BornAnalyzer_cxx
