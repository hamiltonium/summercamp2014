//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 17 00:15:22 2014 by ROOT version 5.32/00
// from TChain mfvFlatTree/t/
//////////////////////////////////////////////////////////

#ifndef MFVFlatNtupleReader_h
#define MFVFlatNtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class MFVFlatNtupleReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;
   Char_t          sample;
   Bool_t          gen_valid;
   Float_t         gen_lsp_pt[2];
   Float_t         gen_lsp_eta[2];
   Float_t         gen_lsp_phi[2];
   Float_t         gen_lsp_mass[2];
   Float_t         gen_lsp_decay[6];
   UChar_t         gen_decay_type[2];
   UChar_t         gen_partons_in_acc;
   Float_t         npu;
   Float_t         bsx;
   Float_t         bsy;
   Float_t         bsz;
   UChar_t         npv;
   Float_t         pvx;
   Float_t         pvy;
   Float_t         pvz;
   UChar_t         pv_ntracks;
   Float_t         pv_sumpt2;
   std::vector<unsigned char> *jet_id;
   std::vector<float>   *jet_pt;
   std::vector<float>   *jet_eta;
   std::vector<float>   *jet_phi;
   std::vector<float>   *jet_energy;
   Float_t         metx;
   Float_t         mety;
   Float_t         metsig;
   Float_t         metdphimin;
   std::vector<unsigned char> *lep_id;
   std::vector<float>   *lep_pt;
   std::vector<float>   *lep_eta;
   std::vector<float>   *lep_phi;
   std::vector<float>   *lep_dxy;
   std::vector<float>   *lep_dz;
   std::vector<float>   *lep_iso;
   std::vector<float>   *lep_mva;
   UChar_t         nvertices;
   std::vector<unsigned char> *vtx_nmu;
   std::vector<unsigned char> *vtx_nel;
   std::vector<float>   *vtx_x;
   std::vector<float>   *vtx_y;
   std::vector<float>   *vtx_z;
   std::vector<float>   *vtx_cxx;
   std::vector<float>   *vtx_cxy;
   std::vector<float>   *vtx_cxz;
   std::vector<float>   *vtx_cyy;
   std::vector<float>   *vtx_cyz;
   std::vector<float>   *vtx_czz;
   std::vector<float>   *vtx_chi2;
   std::vector<float>   *vtx_ndof;
   std::vector<unsigned char> *vtx_njets;
   std::vector<float>   *vtx_tks_pt;
   std::vector<float>   *vtx_tks_eta;
   std::vector<float>   *vtx_tks_phi;
   std::vector<float>   *vtx_tks_mass;
   std::vector<float>   *vtx_jets_pt;
   std::vector<float>   *vtx_jets_eta;
   std::vector<float>   *vtx_jets_phi;
   std::vector<float>   *vtx_jets_mass;
   std::vector<float>   *vtx_tksjets_pt;
   std::vector<float>   *vtx_tksjets_eta;
   std::vector<float>   *vtx_tksjets_phi;
   std::vector<float>   *vtx_tksjets_mass;
   std::vector<float>   *vtx_jetpairdetamin;
   std::vector<float>   *vtx_jetpairdetamax;
   std::vector<float>   *vtx_jetpairdetaavg;
   std::vector<float>   *vtx_jetpairdetarms;
   std::vector<float>   *vtx_jetpairdrmin;
   std::vector<float>   *vtx_jetpairdrmax;
   std::vector<float>   *vtx_jetpairdravg;
   std::vector<float>   *vtx_jetpairdrrms;
   std::vector<float>   *vtx_costhtkmomvtxdispmin;
   std::vector<float>   *vtx_costhtkmomvtxdispmax;
   std::vector<float>   *vtx_costhtkmomvtxdispavg;
   std::vector<float>   *vtx_costhtkmomvtxdisprms;
   std::vector<float>   *vtx_costhjetmomvtxdispmin;
   std::vector<float>   *vtx_costhjetmomvtxdispmax;
   std::vector<float>   *vtx_costhjetmomvtxdispavg;
   std::vector<float>   *vtx_costhjetmomvtxdisprms;
   std::vector<float>   *vtx_gen2ddist;
   std::vector<float>   *vtx_gen2derr;
   std::vector<float>   *vtx_gen3ddist;
   std::vector<float>   *vtx_gen3derr;
   std::vector<float>   *vtx_bs2ddist;
   std::vector<float>   *vtx_bs2derr;
   std::vector<float>   *vtx_pv2ddist;
   std::vector<float>   *vtx_pv2derr;
   std::vector<float>   *vtx_pv3ddist;
   std::vector<float>   *vtx_pv3derr;
   std::vector<unsigned char> *vtx_ntracks;
   std::vector<unsigned char> *vtx_nbadtracks;
   std::vector<unsigned char> *vtx_ntracksptgt3;
   std::vector<unsigned char> *vtx_ntracksptgt5;
   std::vector<unsigned char> *vtx_ntracksptgt10;
   std::vector<unsigned char> *vtx_trackminnhits;
   std::vector<unsigned char> *vtx_trackmaxnhits;
   std::vector<float>   *vtx_sumpt2;
   std::vector<unsigned char> *vtx_sumnhitsbehind;
   std::vector<unsigned char> *vtx_maxnhitsbehind;
   std::vector<unsigned char> *vtx_ntrackssharedwpv;
   std::vector<unsigned char> *vtx_ntrackssharedwpvs;
   std::vector<unsigned char> *vtx_npvswtracksshared;
   std::vector<char>    *vtx_pvmosttracksshared;
   std::vector<float>   *vtx_mintrackpt;
   std::vector<float>   *vtx_maxtrackpt;
   std::vector<float>   *vtx_maxm1trackpt;
   std::vector<float>   *vtx_maxm2trackpt;
   std::vector<float>   *vtx_trackpairdrmin;
   std::vector<float>   *vtx_trackpairdrmax;
   std::vector<float>   *vtx_trackpairdravg;
   std::vector<float>   *vtx_trackpairdrrms;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_sample;   //!
   TBranch        *b_gen_valid;   //!
   TBranch        *b_gen_lsp_pt;   //!
   TBranch        *b_gen_lsp_eta;   //!
   TBranch        *b_gen_lsp_phi;   //!
   TBranch        *b_gen_lsp_mass;   //!
   TBranch        *b_gen_lsp_decay;   //!
   TBranch        *b_gen_decay_type;   //!
   TBranch        *b_gen_partons_in_acc;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_bsx;   //!
   TBranch        *b_bsy;   //!
   TBranch        *b_bsz;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_pvx;   //!
   TBranch        *b_pvy;   //!
   TBranch        *b_pvz;   //!
   TBranch        *b_pv_ntracks;   //!
   TBranch        *b_pv_sumpt2;   //!
   TBranch        *b_jet_id;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_metx;   //!
   TBranch        *b_mety;   //!
   TBranch        *b_metsig;   //!
   TBranch        *b_metdphimin;   //!
   TBranch        *b_lep_id;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_dxy;   //!
   TBranch        *b_lep_dz;   //!
   TBranch        *b_lep_iso;   //!
   TBranch        *b_lep_mva;   //!
   TBranch        *b_nvertices;   //!
   TBranch        *b_vtx_nmu;   //!
   TBranch        *b_vtx_nel;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_cxx;   //!
   TBranch        *b_vtx_cxy;   //!
   TBranch        *b_vtx_cxz;   //!
   TBranch        *b_vtx_cyy;   //!
   TBranch        *b_vtx_cyz;   //!
   TBranch        *b_vtx_czz;   //!
   TBranch        *b_vtx_chi2;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_njets;   //!
   TBranch        *b_vtx_tks_pt;   //!
   TBranch        *b_vtx_tks_eta;   //!
   TBranch        *b_vtx_tks_phi;   //!
   TBranch        *b_vtx_tks_mass;   //!
   TBranch        *b_vtx_jets_pt;   //!
   TBranch        *b_vtx_jets_eta;   //!
   TBranch        *b_vtx_jets_phi;   //!
   TBranch        *b_vtx_jets_mass;   //!
   TBranch        *b_vtx_tksjets_pt;   //!
   TBranch        *b_vtx_tksjets_eta;   //!
   TBranch        *b_vtx_tksjets_phi;   //!
   TBranch        *b_vtx_tksjets_mass;   //!
   TBranch        *b_vtx_jetpairdetamin;   //!
   TBranch        *b_vtx_jetpairdetamax;   //!
   TBranch        *b_vtx_jetpairdetaavg;   //!
   TBranch        *b_vtx_jetpairdetarms;   //!
   TBranch        *b_vtx_jetpairdrmin;   //!
   TBranch        *b_vtx_jetpairdrmax;   //!
   TBranch        *b_vtx_jetpairdravg;   //!
   TBranch        *b_vtx_jetpairdrrms;   //!
   TBranch        *b_vtx_costhtkmomvtxdispmin;   //!
   TBranch        *b_vtx_costhtkmomvtxdispmax;   //!
   TBranch        *b_vtx_costhtkmomvtxdispavg;   //!
   TBranch        *b_vtx_costhtkmomvtxdisprms;   //!
   TBranch        *b_vtx_costhjetmomvtxdispmin;   //!
   TBranch        *b_vtx_costhjetmomvtxdispmax;   //!
   TBranch        *b_vtx_costhjetmomvtxdispavg;   //!
   TBranch        *b_vtx_costhjetmomvtxdisprms;   //!
   TBranch        *b_vtx_gen2ddist;   //!
   TBranch        *b_vtx_gen2derr;   //!
   TBranch        *b_vtx_gen3ddist;   //!
   TBranch        *b_vtx_gen3derr;   //!
   TBranch        *b_vtx_bs2ddist;   //!
   TBranch        *b_vtx_bs2derr;   //!
   TBranch        *b_vtx_pv2ddist;   //!
   TBranch        *b_vtx_pv2derr;   //!
   TBranch        *b_vtx_pv3ddist;   //!
   TBranch        *b_vtx_pv3derr;   //!
   TBranch        *b_vtx_ntracks;   //!
   TBranch        *b_vtx_nbadtracks;   //!
   TBranch        *b_vtx_ntracksptgt3;   //!
   TBranch        *b_vtx_ntracksptgt5;   //!
   TBranch        *b_vtx_ntracksptgt10;   //!
   TBranch        *b_vtx_trackminnhits;   //!
   TBranch        *b_vtx_trackmaxnhits;   //!
   TBranch        *b_vtx_sumpt2;   //!
   TBranch        *b_vtx_sumnhitsbehind;   //!
   TBranch        *b_vtx_maxnhitsbehind;   //!
   TBranch        *b_vtx_ntrackssharedwpv;   //!
   TBranch        *b_vtx_ntrackssharedwpvs;   //!
   TBranch        *b_vtx_npvswtracksshared;   //!
   TBranch        *b_vtx_pvmosttracksshared;   //!
   TBranch        *b_vtx_mintrackpt;   //!
   TBranch        *b_vtx_maxtrackpt;   //!
   TBranch        *b_vtx_maxm1trackpt;   //!
   TBranch        *b_vtx_maxm2trackpt;   //!
   TBranch        *b_vtx_trackpairdrmin;   //!
   TBranch        *b_vtx_trackpairdrmax;   //!
   TBranch        *b_vtx_trackpairdravg;   //!
   TBranch        *b_vtx_trackpairdrrms;   //!

   MFVFlatNtupleReader(TTree *tree=0);
   virtual ~MFVFlatNtupleReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MFVFlatNtupleReader_cxx
MFVFlatNtupleReader::MFVFlatNtupleReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("mfvFlatTree/t",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("mfvFlatTree/t","");
#if defined(ATLPC)
      TString path = "root://cmseos.fnal.gov//store/user/tucker/flattrees/";
#elif defined(ATLNX)
      TString path = "/cdat/tem/jmt46/flattrees/";
#else
      TString path = "";
#endif
      chain->Add(path + "mfv_neutralino_tau0100um_M0400.root/mfvFlatTree/t");
      chain->Add(path + "mfv_neutralino_tau0300um_M0400.root/mfvFlatTree/t");
      chain->Add(path + "mfv_neutralino_tau1000um_M0400.root/mfvFlatTree/t");
      chain->Add(path + "mfv_neutralino_tau9900um_M0400.root/mfvFlatTree/t");
      chain->Add(path + "qcdht0100.root/mfvFlatTree/t");
      chain->Add(path + "qcdht0250.root/mfvFlatTree/t");
      chain->Add(path + "qcdht0500.root/mfvFlatTree/t");
      chain->Add(path + "qcdht1000.root/mfvFlatTree/t");
      chain->Add(path + "ttbardilep.root/mfvFlatTree/t");
      chain->Add(path + "ttbarhadronic.root/mfvFlatTree/t");
      chain->Add(path + "ttbarsemilep.root/mfvFlatTree/t");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

MFVFlatNtupleReader::~MFVFlatNtupleReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MFVFlatNtupleReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MFVFlatNtupleReader::LoadTree(Long64_t entry)
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

void MFVFlatNtupleReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jet_id = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_energy = 0;
   lep_id = 0;
   lep_pt = 0;
   lep_eta = 0;
   lep_phi = 0;
   lep_dxy = 0;
   lep_dz = 0;
   lep_iso = 0;
   lep_mva = 0;
   vtx_nmu = 0;
   vtx_nel = 0;
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   vtx_cxx = 0;
   vtx_cxy = 0;
   vtx_cxz = 0;
   vtx_cyy = 0;
   vtx_cyz = 0;
   vtx_czz = 0;
   vtx_chi2 = 0;
   vtx_ndof = 0;
   vtx_njets = 0;
   vtx_tks_pt = 0;
   vtx_tks_eta = 0;
   vtx_tks_phi = 0;
   vtx_tks_mass = 0;
   vtx_jets_pt = 0;
   vtx_jets_eta = 0;
   vtx_jets_phi = 0;
   vtx_jets_mass = 0;
   vtx_tksjets_pt = 0;
   vtx_tksjets_eta = 0;
   vtx_tksjets_phi = 0;
   vtx_tksjets_mass = 0;
   vtx_jetpairdetamin = 0;
   vtx_jetpairdetamax = 0;
   vtx_jetpairdetaavg = 0;
   vtx_jetpairdetarms = 0;
   vtx_jetpairdrmin = 0;
   vtx_jetpairdrmax = 0;
   vtx_jetpairdravg = 0;
   vtx_jetpairdrrms = 0;
   vtx_costhtkmomvtxdispmin = 0;
   vtx_costhtkmomvtxdispmax = 0;
   vtx_costhtkmomvtxdispavg = 0;
   vtx_costhtkmomvtxdisprms = 0;
   vtx_costhjetmomvtxdispmin = 0;
   vtx_costhjetmomvtxdispmax = 0;
   vtx_costhjetmomvtxdispavg = 0;
   vtx_costhjetmomvtxdisprms = 0;
   vtx_gen2ddist = 0;
   vtx_gen2derr = 0;
   vtx_gen3ddist = 0;
   vtx_gen3derr = 0;
   vtx_bs2ddist = 0;
   vtx_bs2derr = 0;
   vtx_pv2ddist = 0;
   vtx_pv2derr = 0;
   vtx_pv3ddist = 0;
   vtx_pv3derr = 0;
   vtx_ntracks = 0;
   vtx_nbadtracks = 0;
   vtx_ntracksptgt3 = 0;
   vtx_ntracksptgt5 = 0;
   vtx_ntracksptgt10 = 0;
   vtx_trackminnhits = 0;
   vtx_trackmaxnhits = 0;
   vtx_sumpt2 = 0;
   vtx_sumnhitsbehind = 0;
   vtx_maxnhitsbehind = 0;
   vtx_ntrackssharedwpv = 0;
   vtx_ntrackssharedwpvs = 0;
   vtx_npvswtracksshared = 0;
   vtx_pvmosttracksshared = 0;
   vtx_mintrackpt = 0;
   vtx_maxtrackpt = 0;
   vtx_maxm1trackpt = 0;
   vtx_maxm2trackpt = 0;
   vtx_trackpairdrmin = 0;
   vtx_trackpairdrmax = 0;
   vtx_trackpairdravg = 0;
   vtx_trackpairdrrms = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("sample", &sample, &b_sample);
   fChain->SetBranchAddress("gen_valid", &gen_valid, &b_gen_valid);
   fChain->SetBranchAddress("gen_lsp_pt", gen_lsp_pt, &b_gen_lsp_pt);
   fChain->SetBranchAddress("gen_lsp_eta", gen_lsp_eta, &b_gen_lsp_eta);
   fChain->SetBranchAddress("gen_lsp_phi", gen_lsp_phi, &b_gen_lsp_phi);
   fChain->SetBranchAddress("gen_lsp_mass", gen_lsp_mass, &b_gen_lsp_mass);
   fChain->SetBranchAddress("gen_lsp_decay", gen_lsp_decay, &b_gen_lsp_decay);
   fChain->SetBranchAddress("gen_decay_type", gen_decay_type, &b_gen_decay_type);
   fChain->SetBranchAddress("gen_partons_in_acc", &gen_partons_in_acc, &b_gen_partons_in_acc);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("bsx", &bsx, &b_bsx);
   fChain->SetBranchAddress("bsy", &bsy, &b_bsy);
   fChain->SetBranchAddress("bsz", &bsz, &b_bsz);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("pvx", &pvx, &b_pvx);
   fChain->SetBranchAddress("pvy", &pvy, &b_pvy);
   fChain->SetBranchAddress("pvz", &pvz, &b_pvz);
   fChain->SetBranchAddress("pv_ntracks", &pv_ntracks, &b_pv_ntracks);
   fChain->SetBranchAddress("pv_sumpt2", &pv_sumpt2, &b_pv_sumpt2);
   fChain->SetBranchAddress("jet_id", &jet_id, &b_jet_id);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("metx", &metx, &b_metx);
   fChain->SetBranchAddress("mety", &mety, &b_mety);
   fChain->SetBranchAddress("metsig", &metsig, &b_metsig);
   fChain->SetBranchAddress("metdphimin", &metdphimin, &b_metdphimin);
   fChain->SetBranchAddress("lep_id", &lep_id, &b_lep_id);
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_dxy", &lep_dxy, &b_lep_dxy);
   fChain->SetBranchAddress("lep_dz", &lep_dz, &b_lep_dz);
   fChain->SetBranchAddress("lep_iso", &lep_iso, &b_lep_iso);
   fChain->SetBranchAddress("lep_mva", &lep_mva, &b_lep_mva);
   fChain->SetBranchAddress("nvertices", &nvertices, &b_nvertices);
   fChain->SetBranchAddress("vtx_nmu", &vtx_nmu, &b_vtx_nmu);
   fChain->SetBranchAddress("vtx_nel", &vtx_nel, &b_vtx_nel);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_cxx", &vtx_cxx, &b_vtx_cxx);
   fChain->SetBranchAddress("vtx_cxy", &vtx_cxy, &b_vtx_cxy);
   fChain->SetBranchAddress("vtx_cxz", &vtx_cxz, &b_vtx_cxz);
   fChain->SetBranchAddress("vtx_cyy", &vtx_cyy, &b_vtx_cyy);
   fChain->SetBranchAddress("vtx_cyz", &vtx_cyz, &b_vtx_cyz);
   fChain->SetBranchAddress("vtx_czz", &vtx_czz, &b_vtx_czz);
   fChain->SetBranchAddress("vtx_chi2", &vtx_chi2, &b_vtx_chi2);
   fChain->SetBranchAddress("vtx_ndof", &vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_njets", &vtx_njets, &b_vtx_njets);
   fChain->SetBranchAddress("vtx_tks_pt", &vtx_tks_pt, &b_vtx_tks_pt);
   fChain->SetBranchAddress("vtx_tks_eta", &vtx_tks_eta, &b_vtx_tks_eta);
   fChain->SetBranchAddress("vtx_tks_phi", &vtx_tks_phi, &b_vtx_tks_phi);
   fChain->SetBranchAddress("vtx_tks_mass", &vtx_tks_mass, &b_vtx_tks_mass);
   fChain->SetBranchAddress("vtx_jets_pt", &vtx_jets_pt, &b_vtx_jets_pt);
   fChain->SetBranchAddress("vtx_jets_eta", &vtx_jets_eta, &b_vtx_jets_eta);
   fChain->SetBranchAddress("vtx_jets_phi", &vtx_jets_phi, &b_vtx_jets_phi);
   fChain->SetBranchAddress("vtx_jets_mass", &vtx_jets_mass, &b_vtx_jets_mass);
   fChain->SetBranchAddress("vtx_tksjets_pt", &vtx_tksjets_pt, &b_vtx_tksjets_pt);
   fChain->SetBranchAddress("vtx_tksjets_eta", &vtx_tksjets_eta, &b_vtx_tksjets_eta);
   fChain->SetBranchAddress("vtx_tksjets_phi", &vtx_tksjets_phi, &b_vtx_tksjets_phi);
   fChain->SetBranchAddress("vtx_tksjets_mass", &vtx_tksjets_mass, &b_vtx_tksjets_mass);
   fChain->SetBranchAddress("vtx_jetpairdetamin", &vtx_jetpairdetamin, &b_vtx_jetpairdetamin);
   fChain->SetBranchAddress("vtx_jetpairdetamax", &vtx_jetpairdetamax, &b_vtx_jetpairdetamax);
   fChain->SetBranchAddress("vtx_jetpairdetaavg", &vtx_jetpairdetaavg, &b_vtx_jetpairdetaavg);
   fChain->SetBranchAddress("vtx_jetpairdetarms", &vtx_jetpairdetarms, &b_vtx_jetpairdetarms);
   fChain->SetBranchAddress("vtx_jetpairdrmin", &vtx_jetpairdrmin, &b_vtx_jetpairdrmin);
   fChain->SetBranchAddress("vtx_jetpairdrmax", &vtx_jetpairdrmax, &b_vtx_jetpairdrmax);
   fChain->SetBranchAddress("vtx_jetpairdravg", &vtx_jetpairdravg, &b_vtx_jetpairdravg);
   fChain->SetBranchAddress("vtx_jetpairdrrms", &vtx_jetpairdrrms, &b_vtx_jetpairdrrms);
   fChain->SetBranchAddress("vtx_costhtkmomvtxdispmin", &vtx_costhtkmomvtxdispmin, &b_vtx_costhtkmomvtxdispmin);
   fChain->SetBranchAddress("vtx_costhtkmomvtxdispmax", &vtx_costhtkmomvtxdispmax, &b_vtx_costhtkmomvtxdispmax);
   fChain->SetBranchAddress("vtx_costhtkmomvtxdispavg", &vtx_costhtkmomvtxdispavg, &b_vtx_costhtkmomvtxdispavg);
   fChain->SetBranchAddress("vtx_costhtkmomvtxdisprms", &vtx_costhtkmomvtxdisprms, &b_vtx_costhtkmomvtxdisprms);
   fChain->SetBranchAddress("vtx_costhjetmomvtxdispmin", &vtx_costhjetmomvtxdispmin, &b_vtx_costhjetmomvtxdispmin);
   fChain->SetBranchAddress("vtx_costhjetmomvtxdispmax", &vtx_costhjetmomvtxdispmax, &b_vtx_costhjetmomvtxdispmax);
   fChain->SetBranchAddress("vtx_costhjetmomvtxdispavg", &vtx_costhjetmomvtxdispavg, &b_vtx_costhjetmomvtxdispavg);
   fChain->SetBranchAddress("vtx_costhjetmomvtxdisprms", &vtx_costhjetmomvtxdisprms, &b_vtx_costhjetmomvtxdisprms);
   fChain->SetBranchAddress("vtx_gen2ddist", &vtx_gen2ddist, &b_vtx_gen2ddist);
   fChain->SetBranchAddress("vtx_gen2derr", &vtx_gen2derr, &b_vtx_gen2derr);
   fChain->SetBranchAddress("vtx_gen3ddist", &vtx_gen3ddist, &b_vtx_gen3ddist);
   fChain->SetBranchAddress("vtx_gen3derr", &vtx_gen3derr, &b_vtx_gen3derr);
   fChain->SetBranchAddress("vtx_bs2ddist", &vtx_bs2ddist, &b_vtx_bs2ddist);
   fChain->SetBranchAddress("vtx_bs2derr", &vtx_bs2derr, &b_vtx_bs2derr);
   fChain->SetBranchAddress("vtx_pv2ddist", &vtx_pv2ddist, &b_vtx_pv2ddist);
   fChain->SetBranchAddress("vtx_pv2derr", &vtx_pv2derr, &b_vtx_pv2derr);
   fChain->SetBranchAddress("vtx_pv3ddist", &vtx_pv3ddist, &b_vtx_pv3ddist);
   fChain->SetBranchAddress("vtx_pv3derr", &vtx_pv3derr, &b_vtx_pv3derr);
   fChain->SetBranchAddress("vtx_ntracks", &vtx_ntracks, &b_vtx_ntracks);
   fChain->SetBranchAddress("vtx_nbadtracks", &vtx_nbadtracks, &b_vtx_nbadtracks);
   fChain->SetBranchAddress("vtx_ntracksptgt3", &vtx_ntracksptgt3, &b_vtx_ntracksptgt3);
   fChain->SetBranchAddress("vtx_ntracksptgt5", &vtx_ntracksptgt5, &b_vtx_ntracksptgt5);
   fChain->SetBranchAddress("vtx_ntracksptgt10", &vtx_ntracksptgt10, &b_vtx_ntracksptgt10);
   fChain->SetBranchAddress("vtx_trackminnhits", &vtx_trackminnhits, &b_vtx_trackminnhits);
   fChain->SetBranchAddress("vtx_trackmaxnhits", &vtx_trackmaxnhits, &b_vtx_trackmaxnhits);
   fChain->SetBranchAddress("vtx_sumpt2", &vtx_sumpt2, &b_vtx_sumpt2);
   fChain->SetBranchAddress("vtx_sumnhitsbehind", &vtx_sumnhitsbehind, &b_vtx_sumnhitsbehind);
   fChain->SetBranchAddress("vtx_maxnhitsbehind", &vtx_maxnhitsbehind, &b_vtx_maxnhitsbehind);
   fChain->SetBranchAddress("vtx_ntrackssharedwpv", &vtx_ntrackssharedwpv, &b_vtx_ntrackssharedwpv);
   fChain->SetBranchAddress("vtx_ntrackssharedwpvs", &vtx_ntrackssharedwpvs, &b_vtx_ntrackssharedwpvs);
   fChain->SetBranchAddress("vtx_npvswtracksshared", &vtx_npvswtracksshared, &b_vtx_npvswtracksshared);
   fChain->SetBranchAddress("vtx_pvmosttracksshared", &vtx_pvmosttracksshared, &b_vtx_pvmosttracksshared);
   fChain->SetBranchAddress("vtx_mintrackpt", &vtx_mintrackpt, &b_vtx_mintrackpt);
   fChain->SetBranchAddress("vtx_maxtrackpt", &vtx_maxtrackpt, &b_vtx_maxtrackpt);
   fChain->SetBranchAddress("vtx_maxm1trackpt", &vtx_maxm1trackpt, &b_vtx_maxm1trackpt);
   fChain->SetBranchAddress("vtx_maxm2trackpt", &vtx_maxm2trackpt, &b_vtx_maxm2trackpt);
   fChain->SetBranchAddress("vtx_trackpairdrmin", &vtx_trackpairdrmin, &b_vtx_trackpairdrmin);
   fChain->SetBranchAddress("vtx_trackpairdrmax", &vtx_trackpairdrmax, &b_vtx_trackpairdrmax);
   fChain->SetBranchAddress("vtx_trackpairdravg", &vtx_trackpairdravg, &b_vtx_trackpairdravg);
   fChain->SetBranchAddress("vtx_trackpairdrrms", &vtx_trackpairdrrms, &b_vtx_trackpairdrrms);
   Notify();
}

Bool_t MFVFlatNtupleReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MFVFlatNtupleReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MFVFlatNtupleReader::Cut(Long64_t)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MFVFlatNtupleReader_cxx
