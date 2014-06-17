#define MFVFlatNtupleReader_cxx
#include "MFVFlatNtupleReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>

void MFVFlatNtupleReader::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L MFVFlatNtupleReader.C
  //      Root > MFVFlatNtupleReader t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
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
  if (fChain == 0) return;


  enum {
    s_start = -3,
    s_mfv_neutralino_tau0100um_M0400 = -3,
    s_mfv_neutralino_tau0300um_M0400,
    s_mfv_neutralino_tau1000um_M0400,
    s_mfv_neutralino_tau9900um_M0400,
    s_ttbarhadronic,
    s_ttbarsemilep,
    s_ttbardilep,
    s_qcdht0100,
    s_qcdht0250,
    s_qcdht0500,
    s_qcdht1000,
    s_end
  };

  std::map<int, std::string> sample_names;
  sample_names[s_mfv_neutralino_tau0100um_M0400] = "mfv_neutralino_tau0100um_M0400";
  sample_names[s_mfv_neutralino_tau0300um_M0400] = "mfv_neutralino_tau0300um_M0400";
  sample_names[s_mfv_neutralino_tau1000um_M0400] = "mfv_neutralino_tau1000um_M0400";
  sample_names[s_mfv_neutralino_tau9900um_M0400] = "mfv_neutralino_tau9900um_M0400";
  sample_names[s_ttbarhadronic] = "ttbarhadronic";
  sample_names[s_ttbarsemilep] = "ttbarsemilep";
  sample_names[s_ttbardilep] = "ttbardilep";
  sample_names[s_qcdht0100] = "qcdht0100";
  sample_names[s_qcdht0250] = "qcdht0250";
  sample_names[s_qcdht0500] = "qcdht0500";
  sample_names[s_qcdht1000] = "qcdht1000";

  std::map<int, int> sample_nevents_read;
  std::map<int, int> sample_nevents_1vtx;
  std::map<int, int> sample_nevents_2vtx;
  std::map<int, int> sample_nevents_2vtx500um;

  Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    if (jentry % 1000 == 0) {
      printf("\r%lli/%lli events read", jentry, nentries);
      fflush(stdout);
    }

    // if (Cut(ientry) < 0) continue;

    ++sample_nevents_read[sample];

    std::vector<size_t> vtxpass;
    for (size_t ivtx = 0; ivtx < nvertices; ++ivtx) {
      if (vtx_ntracks->at(ivtx) >= 5 &&
          vtx_trackpairdrmin->at(ivtx) < 0.4 &&
          vtx_trackpairdrmax->at(ivtx) > 1.2 &&
          vtx_trackpairdrmax->at(ivtx) > 4 &&
          vtx_bs2derr->at(ivtx) < 0.0025 &&
          vtx_njets->at(ivtx) >= 1 &&
          vtx_ntracksptgt3->at(ivtx) >= 3 &&
          vtx_sumnhitsbehind->at(ivtx) == 0) {
        vtxpass.push_back(ivtx);
      }
    }

    if (vtxpass.size() >= 1) ++sample_nevents_1vtx[sample];
    if (vtxpass.size() >= 2) {
      ++sample_nevents_2vtx[sample];

      bool gt500um = false;
      for (size_t ii = 0; ii < vtxpass.size(); ++ii) {
        size_t ivtx = vtxpass[ii];
        for (size_t jj = ii+1; jj < vtxpass.size(); ++jj) {
          size_t jvtx = vtxpass[jj];
          if (sqrt(pow(vtx_x->at(ivtx) - vtx_x->at(jvtx), 2) +
                   pow(vtx_y->at(ivtx) - vtx_y->at(jvtx), 2) +
                   pow(vtx_z->at(ivtx) - vtx_z->at(jvtx), 2)) > 0.05)
            gt500um = true;
        }
      }

      if (gt500um)
        ++sample_nevents_2vtx500um[sample];
    }
  }

  printf("\r%80s\rdone!\n", "");

  for (int s = s_start; s < s_end; ++s)
    printf("sample %30s events read: %10i\n", sample_names[s].c_str(), sample_nevents_read[s]);

  printf("\n");
  printf("qcdht1000    events read: %10i\n", sample_nevents_read[s_qcdht1000]);
  printf("             # 1-vtx evs: %10i\n", sample_nevents_1vtx[s_qcdht1000]);
  printf("             # 2-vtx evs: %10i\n", sample_nevents_2vtx[s_qcdht1000]);
  printf("# 2-vtx sep by 500um evs: %10i\n", sample_nevents_2vtx500um[s_qcdht1000]);
}

#ifdef STANDALONE
int main() {
  MFVFlatNtupleReader r;
  r.Loop();
}
#endif
