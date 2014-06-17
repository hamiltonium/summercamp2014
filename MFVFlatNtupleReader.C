#define MFVFlatNtupleReader_cxx
#include "MFVFlatNtupleReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
  }

  printf("\r%80s\rdone!\n", "");

  for (int s = s_start; s < s_end; ++s)
    printf("sample %30s events read: %10i\n", sample_names[s].c_str(), sample_nevents_read[s]);
}

#ifdef STANDALONE
int main() {
  MFVFlatNtupleReader r;
  r.Loop();
}
#endif
