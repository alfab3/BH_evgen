#define BornAnalyzer_cxx
#include "BornAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

using namespace ROOT::Math;
const double ElectronMass = 5.109989461e-4; //Mass of electron


//put some BS here lol

XYZVector ThreeVec(double Theta, double Phi)
{
  XYZVector k(  ( sin(Theta) * cos(Phi)),  sin(Theta)*sin(Phi), cos(Theta) ) ; //three momentum for lepton
  return k;
}


void BornAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L BornAnalyzer.C
//      root> BornAnalyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   TLorentzVector p1; TLorentzVector p2; TLorentzVector BeamP4;

   TFile* oFile = TFile::Open("Born2BHHists.root", "RECREATE");
   TH1D* Hist_W2e = new TH1D("W2e", ";Invariant Mass (GeV/c^2)", 234, 0, 1.4);
   TH1D* Hist_t = new TH1D("t", ";Momentum Transfer Squared", 234, 0, .14 );


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;



 XYZVector k1 = ThreeVec(theta1, phi1) * sqrt( pow(x,2) * pow(E0,2) - pow(ElectronMass,2)); 
      XYZVector k2 = ThreeVec(theta2, phi2) * sqrt( pow(1 - x,2) * pow(E0,2) - pow(ElectronMass,2));

      p1.SetXYZM(k1.X(),k1.Y(),k1.Z(),ElectronMass); //v=(x,y,z,e=Sqrt(x*x+y*y+z*z+m*m))
      p2.SetXYZM(k2.X(),k2.Y(),k2.Z(),ElectronMass);

      BeamP4.SetXYZT(0,0,E0,E0);

      double W2e = sqrt((p1 + p2).Mag2()); 
      Hist_W2e->Fill(W2e);   

      double t = (BeamP4 - p1 - p2).Mag2();
      Hist_t->Fill(-t);
		//the momentum transfers are really transverse momentum transfers. And you want to have a consistency
		//between how the momentum transfers are handled in the form factors
		//the singularity at forward angles will go away when you put in the atomic form factors



   }
oFile->Write(); 
oFile->Close();
}
