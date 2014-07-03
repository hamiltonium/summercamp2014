#define MFVFlatNtupleReader_cxx
#include "MFVFlatNtupleReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <THStack.h>
#include <TLegend.h>
#include <string.h>

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
    s_mfv_neutralino_tau1000um_M0400,
    s_mfv_neutralino_tau0300um_M0400,
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
    //std::cout << std::endl << b_vtx_sumnhitsbehind->GetEntry(jentry) << std::endl;
    
    int numverts = nvertices;
    int numgoodverts = 0;
        //std::cout << "numverts = " << numverts  << std::endl;
    for (int i = 0; i < numverts; ++i) {
      if (MFVFlatNtupleReader::goodVertex(i)) {
	++(numgoodverts);
	//std::cout << "good vertex detected" << std::endl;
      }
    }
    
    //std::cout << numgoodverts << std::endl;
    if (numgoodverts >= 1) {
      ++sample_nevents_read[sample];
      if (sample < 0) {
	//do things here (signal)
      }
      else {
	//do similar but related things here (background)
      }
      //do other related thing regardless
    }
  }

  printf("\r%80s\rdone!\n", "");

  for (int s = s_start; s < s_end; ++s)
    printf("sample %30s events read: %10i\n", sample_names[s].c_str(), sample_nevents_read[s]);
}




void MFVFlatNtupleReader::loopDisplaySigBack()
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
    s_mfv_neutralino_tau1000um_M0400,
    s_mfv_neutralino_tau0300um_M0400,
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
  
  //variables to check signal to background before and after cuts
  int ngoodrsigevents = 0;
  int ngoodrbackevents = 0;
  
  int ngoodocsigevents = 0;
  int ngoodocbackevents = 0;

  int ngoodncsigevents = 0;
  int ngoodncbackevents = 0;
  
  Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    if (jentry % 1000 == 0) {
      printf("\r%lli/%lli events read", jentry, nentries);
      fflush(stdout);
    }
    //std::cout << std::endl << b_vtx_sumnhitsbehind->GetEntry(jentry) << std::endl;
    
    int numverts = nvertices;
    int numgoodncverts = 0;
    int numgoodocverts = 0;
    //std::cout << "numverts = " << numverts  << std::endl;
    for (int i = 0; i < numverts; ++i) {
      if (MFVFlatNtupleReader::goodVertex(i)) {
	++(numgoodncverts);
	//std::cout << "good vertex detected" << std::endl;
      }
      if (MFVFlatNtupleReader::goodVertex_old(i)) {
	++(numgoodocverts);
      }
    }
    
    //std::cout << numgoodverts << std::endl;
    if (numgoodncverts >= 1) {
      ++sample_nevents_read[sample];
      if (sample < 0)  ++ngoodncsigevents;
      else ++ngoodncbackevents;
    }
    if (numgoodocverts >= 1) {
      if (sample < 0)  ++ngoodocsigevents;
      else ++ngoodocbackevents;
    }
    if (sample < 0)  ++ngoodrsigevents;
    else ++ngoodrbackevents;
  }

  printf("\r%80s\rdone!\n", "");

  for (int s = s_start; s < s_end; ++s)
    printf("sample %30s events read: %10i\n", sample_names[s].c_str(), sample_nevents_read[s]);

  //signal vs background
  std::cout << std::endl << "new cuts: " << std::endl << "good signal events: " << ngoodncsigevents << std::endl;
  std::cout << "good background events: " << ngoodncbackevents << std::endl << "ratio: " << (float)ngoodncsigevents/(float)ngoodncbackevents;
  std::cout << std::endl << std::endl;

  std::cout << std::endl << "old cuts: " << std::endl << "good signal events: " << ngoodocsigevents << std::endl;
  std::cout << "good background events: " << ngoodocbackevents << std::endl << "ratio: " << (float)ngoodocsigevents/(float)ngoodocbackevents;
  std::cout << std::endl << std::endl;

  std::cout << std::endl << "no cuts: " << std::endl << "good signal events: " << ngoodrsigevents << std::endl;
  std::cout << "good background events: " << ngoodrbackevents << std::endl << "ratio: " << (float)ngoodrsigevents/(float)ngoodrbackevents;
  std::cout << std::endl << std::endl;
}

void MFVFlatNtupleReader::genVtx_ntracksHists()
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
  
  //set up the histograms to desired specifications
  int nbins = 10;
  h_vtx_ntracks_tot = new TH1F("h_vtx_ntracks_tot","Vtx_ntracks After Cut (Total)",nbins,0.5,nbins+0.5);
  h_vtx_ntracks_sig = new TH1F("h_vtx_ntracks_sig","Vtx_ntracks After Cut (Signal)",nbins,0.5,nbins+0.5);
  h_vtx_ntracks_back = new TH1F("h_vtx_ntracks_back","Vtx_ntracks After Cut (Background)",nbins,0.5,nbins+0.5);
  
  enum {
    s_start = -3,
    s_mfv_neutralino_tau0100um_M0400 = -3,
    s_mfv_neutralino_tau1000um_M0400,
    s_mfv_neutralino_tau0300um_M0400,
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
    //std::cout << std::endl << b_vtx_sumnhitsbehind->GetEntry(jentry) << std::endl;
    
    int numverts = nvertices;
    int numgoodverts = 0;
        //std::cout << "numverts = " << numverts  << std::endl;
    for (int i = 0; i < numverts; ++i) {
      if (MFVFlatNtupleReader::goodVertex(i)) {
	++(numgoodverts);
	//std::cout << "good vertex detected" << std::endl;
      }
    }
    
    //std::cout << numgoodverts << std::endl;
    if (numgoodverts >= 1) {
      ++sample_nevents_read[sample];
      for (int i = 0; i < nvertices; ++i) { //fill with everything
	h_vtx_ntracks_tot->Fill((int)(*vtx_ntracks)[i]);
	  }
      if (sample < 0) { //fill with signal
	for (int i = 0; i < nvertices; ++i) {
	  h_vtx_ntracks_sig->Fill((int)(*vtx_ntracks)[i]);
	  }
      }
      else { //fill with background
	for (int i = 0; i < nvertices; ++i) {
	  h_vtx_ntracks_back->Fill((int)(*vtx_ntracks)[i]);
	  }
      }
    }
  }

  printf("\r%80s\rdone!\n", "");

  for (int s = s_start; s < s_end; ++s)
    printf("sample %30s events read: %10i\n", sample_names[s].c_str(), sample_nevents_read[s]);

  //now we set up the canvas and histogram details so we can draw them
  TCanvas *c1 = new TCanvas("c1","Vtx_ntracks of Lepton-Producing Events",900,600);
  c1->SetBorderMode(0);
  c1->SetFillColor(kWhite);
  //c1->SetLogy();

  THStack *hs = new THStack("hs","Vtx_ntracks of Lepton-Producing Events");
  h_vtx_ntracks_tot->SetFillColor(kGreen);
  hs->Add(h_vtx_ntracks_tot);

  h_vtx_ntracks_back->SetFillColor(kBlue);
  hs->Add(h_vtx_ntracks_back);

  h_vtx_ntracks_sig->SetFillColor(kYellow);
  hs->Add(h_vtx_ntracks_sig);

  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle("vtx_ntracks");
  hs->GetYaxis()->SetTitle("number of vertices");
  hs->GetYaxis()->SetTitleOffset(1.4);
  
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetHeader("Type of Data");
  leg->AddEntry(h_vtx_ntracks_tot,"Total");
  leg->AddEntry(h_vtx_ntracks_back,"Background");
  leg->AddEntry(h_vtx_ntracks_sig,"Signal");
  leg->Draw();

  c1->Modified();
  
  //data to be printed to the command line
  vector<int> entries_in_bins_tot;
  vector<int> entries_in_bins_back;
  vector<int> entries_in_bins_sig;
  for (int i = 1; i <= nbins; ++i)  {
    entries_in_bins_tot.push_back(h_vtx_ntracks_tot->GetBinContent(i));
    entries_in_bins_back.push_back(h_vtx_ntracks_back->GetBinContent(i));
    entries_in_bins_sig.push_back(h_vtx_ntracks_sig->GetBinContent(i));
  }
  
  std::cout << std::endl << "Total: " << std::endl;
  for (int i = 0; i < nbins; ++i) {
    std::cout << "bin " << i+1 <<  ": " << entries_in_bins_tot[i] << " entries" << std::endl;
  }

  std::cout << std::endl << "Background: " << std::endl;
  for (int i = 0; i < nbins; ++i) {
    std::cout << "bin " << i+1 << ": " << entries_in_bins_back[i] << " entries" << std::endl;
  }

  std::cout << std::endl << "Signal: " << std::endl;
  for (int i = 0; i < nbins; ++i) {
    std::cout << "bin " << i+1 << ": " << entries_in_bins_sig[i] <<  " entries" << std::endl;
  }
  std::cout << std::endl;
}

void MFVFlatNtupleReader::genHists(std::string target)
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
  
  //set up the program for a specific target
  std::vector<unsigned char> *loc;
  loc = 0;
  if (target == "vtx_ntracks") fChain->SetBranchAddress("vtx_ntracks", &loc, &b_vtx_ntracks);
  else if (target == "vtx_ntracksptgt3") fChain->SetBranchAddress("vtx_ntracksptgt3",&loc,&b_vtx_ntracksptgt3);
  else if (target == "vtx_ntracksptgt5") fChain->SetBranchAddress("vtx_ntracksptgt5",&loc,&b_vtx_ntracksptgt5);
  else if (target == "vtx_ntracksptgt10") fChain->SetBranchAddress("vtx_ntracksptgt10",&loc,&b_vtx_ntracksptgt10);
  else if (target == "vtx_sumnhitsbehind") fChain->SetBranchAddress("vtx_sumnhitsbehind",&loc,&b_vtx_sumnhitsbehind);
  else if (target == "vtx_njets") fChain->SetBranchAddress("vtx_njets",&loc,&b_vtx_njets);
  else if (target == "vtx_nmu") fChain->SetBranchAddress("vtx_nmu",&loc,&b_vtx_nmu);
  else if (target == "vtx_nel") fChain->SetBranchAddress("vtx_nel",&loc,&b_vtx_nel);
  else if (target == "vtx_maxnhitsbehind") fChain->SetBranchAddress("vtx_maxnhitsbehind",&loc,&b_vtx_sumnhitsbehind);

  //set up the histograms to desired specifications
  int nbins = 26;
  h_tot = new TH1F("h_tot"," After Cut (Total)",nbins,-0.5,nbins-0.5);
  h_back = new TH1F("h_back"," After Cut (Background)",nbins,-0.5,nbins-0.5);
  h_sig = new TH1F("h_sig"," After Cut (Signal)",nbins,-0.5,nbins-0.5);
  
  enum {
    s_start = -3,
    s_mfv_neutralino_tau0100um_M0400 = -3,
    s_mfv_neutralino_tau1000um_M0400,
    s_mfv_neutralino_tau0300um_M0400,
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
    //std::cout << std::endl << b_vtx_sumnhitsbehind->GetEntry(jentry) << std::endl;
    
    int numverts = nvertices;
    int numgoodverts = 0;
        //std::cout << "numverts = " << numverts  << std::endl;
    for (int i = 0; i < numverts; ++i) {
      if (MFVFlatNtupleReader::goodVertex(i)) {
	++(numgoodverts);
	//std::cout << "good vertex detected" << std::endl;
      }
    }
    
    //std::cout << numgoodverts << std::endl;
    if (numgoodverts >= 1) {
      ++sample_nevents_read[sample];
      for (int i = 0; i < nvertices; ++i) { //fill with everything
	h_tot->Fill((int)(*loc)[i]);
	  }
      if (sample >= 0) { //fill with background
	for (int i = 0; i < nvertices; ++i) {
	  h_back->Fill((int)(*loc)[i]);
	}                                                       //ERRORS TO BE FIXED HERE: MAP STRINGS AND POINTERS!!!
      }
      else { //fill with signal
	for (int i = 0; i < nvertices; ++i) {
	  h_sig->Fill((int)(*loc)[i]);
	  }
      }
    }
  }

  printf("\r%80s\rdone!\n", "");

  for (int s = s_start; s < s_end; ++s)
    printf("sample %30s events read: %10i\n", sample_names[s].c_str(), sample_nevents_read[s]);

  gStyle->SetCanvasPreferGL(1);
  //now we set up the canvas and histogram details so we can draw them
  string title_str = target + " of Lepton-Producing Events";
  const char * title = new char[sizeof title_str];
  title = title_str.c_str();
  const char * xaxis = target.c_str();

  TCanvas *c1 = new TCanvas("c1",title,900,600);
  c1->SetBorderMode(0);
  c1->SetFillColor(kWhite);
  //c1->SetLogy();
  
  THStack *hs = new THStack("hs",title);
  h_tot->SetLineColor(kGreen);
  h_tot->SetFillStyle(4050);
  h_tot->SetLineWidth(2);
  hs->Add(h_tot);

  h_back->SetLineColor(kBlue);
  h_back->SetFillStyle(4050);
  h_back->SetLineWidth(2);
  hs->Add(h_back);

  h_sig->SetLineColor(kRed);
  h_sig->SetFillStyle(4010);
  h_sig->SetLineWidth(2);
  hs->Add(h_sig);
  
  //h_tot->Draw();
  //h_back->Draw("same");
  //h_sig->Draw("same");
  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle(xaxis);
  hs->GetYaxis()->SetTitle("number of vertices");
  hs->GetYaxis()->SetTitleOffset(1.4);
  
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetHeader("Type of Data");
  leg->AddEntry(h_tot,"Total");
  leg->AddEntry(h_back,"Background");
  leg->AddEntry(h_sig,"Signal");
  leg->Draw();

  c1->Modified();
  
  //data to be printed to the command line
  vector<int> entries_in_bins_tot;
  vector<int> entries_in_bins_back;
  vector<int> entries_in_bins_sig;
  for (int i = 1; i <= nbins; ++i)  {
    entries_in_bins_tot.push_back(h_tot->GetBinContent(i));
    entries_in_bins_back.push_back(h_back->GetBinContent(i));
    entries_in_bins_sig.push_back(h_sig->GetBinContent(i));
  }
  
  std::cout << std::endl << "Total: " << std::endl;
  for (int i = 0; i < nbins; ++i) {
    std::cout << "bin " << i+1 <<  ": " << entries_in_bins_tot[i] << " entries" << std::endl;
  }

  std::cout << std::endl << "Background: " << std::endl;
  for (int i = 0; i < nbins; ++i) {
    std::cout << "bin " << i+1 << ": " << entries_in_bins_back[i] << " entries" << std::endl;
  }

  std::cout << std::endl << "Signal: " << std::endl;
  for (int i = 0; i < nbins; ++i) {
    std::cout << "bin " << i+1 << ": " << entries_in_bins_sig[i] <<  " entries" << std::endl;
  }
  std::cout << std::endl;
}

/*
void genPercentHists()
{
int nbins = 26;
TH1F * h_back_percent = new TH1F("h_back_percent","Percent of Bin Content (Background)",nbins,-0.5,-0.5+nbins);
TH1F * h_sig_percent = new TH1F("h_sig_percent","Percent of Bin Content (Signal)",nbins,-0.5,-0.5+nbins);

vector<int> entries_in_bins_tot;
  vector<int> entries_in_bins_back;
  vector<int> entries_in_bins_sig;
  for (int i = 1; i <= nbins; ++i)  {
    entries_in_bins_tot.push_back(h_tot->GetBinContent(i));
    entries_in_bins_back.push_back(h_back->GetBinContent(i));
    entries_in_bins_sig.push_back(h_sig->GetBinContent(i));
  }
  for (int i = 0; i < nbins; ++i) {
  if (entries_in_bins_tot[i] == 0) {
  
  }
}
}
 */
#ifdef STANDALONE
int main() {
  MFVFlatNtupleReader r;
  r.Loop();
}
#endif
