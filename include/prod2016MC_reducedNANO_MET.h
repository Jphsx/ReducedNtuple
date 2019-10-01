//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 19 16:43:02 2019 by ROOT version 6.12/07
// from TTree Events/Events
// found on file: prod2016MC_NANO_1.root
//////////////////////////////////////////////////////////

#ifndef prod2016MC_reducedNANO_MET_h
#define prod2016MC_reducedNANO_MET_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class prod2016MC_reducedNANO_MET {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_sumEt;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Bool_t          HLT_PFHT300_PFMET100;
   Bool_t          HLT_PFHT300_PFMET110;
   Bool_t          HLT_PFMET170_NotCleaned;
   Bool_t          HLT_PFMET170_NoiseCleaned;
   Bool_t          HLT_PFMET170_HBHECleaned;
   Bool_t          HLT_PFMET170_JetIdCleaned;
   Bool_t          HLT_PFMET170_BeamHaloCleaned;
   Bool_t          HLT_PFMET170_HBHE_BeamHaloCleaned;
   Bool_t          HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned;
   Bool_t          HLT_PFMET90_PFMHT90_IDTight;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned;
   Bool_t          HLT_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_MET100;
   Bool_t          HLT_MET150;
   Bool_t          HLT_MET200;
   Bool_t          HLT_MET250;
   Bool_t          HLT_MET300;
   Bool_t          HLT_MET600;
   Bool_t          HLT_MET700;
   Bool_t          HLT_PFMET300;
   Bool_t          HLT_PFMET400;
   Bool_t          HLT_PFMET500;
   Bool_t          HLT_PFMET600;

   // List of branches
   TBranch        *b_GenMET_phi;   //!
   TBranch        *b_GenMET_pt;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_HLT_PFHT300_PFMET100;   //!
   TBranch        *b_HLT_PFHT300_PFMET110;   //!
   TBranch        *b_HLT_PFMET170_NotCleaned;   //!
   TBranch        *b_HLT_PFMET170_NoiseCleaned;   //!
   TBranch        *b_HLT_PFMET170_HBHECleaned;   //!
   TBranch        *b_HLT_PFMET170_JetIdCleaned;   //!
   TBranch        *b_HLT_PFMET170_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMET170_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMET90_PFMHT90_IDTight;   //!
   TBranch        *b_HLT_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_MET100;   //!
   TBranch        *b_HLT_MET150;   //!
   TBranch        *b_HLT_MET200;   //!
   TBranch        *b_HLT_MET250;   //!
   TBranch        *b_HLT_MET300;   //!
   TBranch        *b_HLT_MET600;   //!
   TBranch        *b_HLT_MET700;   //!
   TBranch        *b_HLT_PFMET300;   //!
   TBranch        *b_HLT_PFMET400;   //!
   TBranch        *b_HLT_PFMET500;   //!
   TBranch        *b_HLT_PFMET600;   //!

   prod2016MC_reducedNANO_MET(TTree *tree=0);
   virtual ~prod2016MC_reducedNANO_MET();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef prod2016MC_reducedNANO_MET_cxx
prod2016MC_reducedNANO_MET::prod2016MC_reducedNANO_MET(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("prod2016MC_NANO_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("prod2016MC_NANO_1.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

prod2016MC_reducedNANO_MET::~prod2016MC_reducedNANO_MET()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t prod2016MC_reducedNANO_MET::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t prod2016MC_reducedNANO_MET::LoadTree(Long64_t entry)
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

void prod2016MC_reducedNANO_MET::Init(TTree *tree)
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
fChain->SetBranchStatus("*", 0);
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("HLT_PFHT300_PFMET100", &HLT_PFHT300_PFMET100, &b_HLT_PFHT300_PFMET100);
   fChain->SetBranchAddress("HLT_PFHT300_PFMET110", &HLT_PFHT300_PFMET110, &b_HLT_PFHT300_PFMET110);
   fChain->SetBranchAddress("HLT_PFMET170_NotCleaned", &HLT_PFMET170_NotCleaned, &b_HLT_PFMET170_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_NoiseCleaned", &HLT_PFMET170_NoiseCleaned, &b_HLT_PFMET170_NoiseCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned", &HLT_PFMET170_HBHECleaned, &b_HLT_PFMET170_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET170_JetIdCleaned", &HLT_PFMET170_JetIdCleaned, &b_HLT_PFMET170_JetIdCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_BeamHaloCleaned", &HLT_PFMET170_BeamHaloCleaned, &b_HLT_PFMET170_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_HBHE_BeamHaloCleaned", &HLT_PFMET170_HBHE_BeamHaloCleaned, &b_HLT_PFMET170_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned", &HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned, &b_HLT_PFMETTypeOne190_HBHE_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMET90_PFMHT90_IDTight", &HLT_PFMET90_PFMHT90_IDTight, &b_HLT_PFMET90_PFMHT90_IDTight);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight", &HLT_PFMET100_PFMHT100_IDTight, &b_HLT_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned", &HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned, &b_HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_MET100", &HLT_MET100, &b_HLT_MET100);
   fChain->SetBranchAddress("HLT_MET150", &HLT_MET150, &b_HLT_MET150);
   fChain->SetBranchAddress("HLT_MET200", &HLT_MET200, &b_HLT_MET200);
   fChain->SetBranchAddress("HLT_MET250", &HLT_MET250, &b_HLT_MET250);
   fChain->SetBranchAddress("HLT_MET300", &HLT_MET300, &b_HLT_MET300);
   fChain->SetBranchAddress("HLT_MET600", &HLT_MET600, &b_HLT_MET600);
   fChain->SetBranchAddress("HLT_MET700", &HLT_MET700, &b_HLT_MET700);
   fChain->SetBranchAddress("HLT_PFMET300", &HLT_PFMET300, &b_HLT_PFMET300);
   fChain->SetBranchAddress("HLT_PFMET400", &HLT_PFMET400, &b_HLT_PFMET400);
   fChain->SetBranchAddress("HLT_PFMET500", &HLT_PFMET500, &b_HLT_PFMET500);
   fChain->SetBranchAddress("HLT_PFMET600", &HLT_PFMET600, &b_HLT_PFMET600);
   Notify();
}

Bool_t prod2016MC_reducedNANO_MET::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void prod2016MC_reducedNANO_MET::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t prod2016MC_reducedNANO_MET::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef prod2016MC_reducedNANO_MET_cxx
