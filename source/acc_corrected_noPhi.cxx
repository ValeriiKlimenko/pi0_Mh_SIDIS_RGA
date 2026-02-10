// make_data_gen_over_mc.cxx
// Compute corrected yield per bin:  DATA * GEN / MC(rec)
// Also save REC/GEN (acceptance-like) into the same output file.

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>   // std::max

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TAxis.h"

struct Yield {
  double v  = 0.0;  // sum of yields
  double e2 = 0.0;  // sum of error^2
};

static inline uint64_t Key(int xq2, int z) {
  return (uint64_t(uint32_t(xq2)) << 32) | uint32_t(z);
}

static std::unordered_map<uint64_t, Yield>
ReadFitTreeToMap(TFile* f, const char* treeName = "h22_fit")
{
  std::unordered_map<uint64_t, Yield> m;

  if (!f || f->IsZombie()) return m;
  TTree* t = nullptr;
  f->GetObject(treeName, t);
  if (!t) {
    std::cerr << "ERROR: tree '" << treeName << "' not found in " << f->GetName() << "\n";
    return m;
  }

  int    xq2bin = 0;
  int    zbin   = 0;  // "z_pt2_phi_bin"
  double nPions = 0;
  double err    = 0;

  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("xq2bin", 1);
  t->SetBranchStatus("z_pt2_phi_bin", 1);
  t->SetBranchStatus("nPions", 1);
  t->SetBranchStatus("errPions", 1);

  t->SetBranchAddress("xq2bin", &xq2bin);
  t->SetBranchAddress("z_pt2_phi_bin", &zbin);
  t->SetBranchAddress("nPions", &nPions);
  t->SetBranchAddress("errPions", &err);

  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    auto& y = m[Key(xq2bin, zbin)];
    y.v  += nPions;
    y.e2 += err * err;
  }
  return m;
}

int make_data_gen_over_mc(const char* data_fit_file,
                          const char* mc_fit_file,
                          const char* gen_file,
                          const char* out_file = "data_times_gen_over_mc.root",
                          const char* gen_hist_name = "h2_binX_vs_z",
                          const char* fit_tree_name = "h22_fit")
{
  // --- Open inputs ---
  TFile fData(data_fit_file, "READ");
  if (fData.IsZombie()) { std::cerr << "ERROR: cannot open " << data_fit_file << "\n"; return 1; }

  TFile fMC(mc_fit_file, "READ");
  if (fMC.IsZombie()) { std::cerr << "ERROR: cannot open " << mc_fit_file << "\n"; return 2; }

  TFile fGen(gen_file, "READ");
  if (fGen.IsZombie()) { std::cerr << "ERROR: cannot open " << gen_file << "\n"; return 3; }

  TH2D* hGen = nullptr;
  fGen.GetObject(gen_hist_name, hGen);
  if (!hGen) {
    std::cerr << "ERROR: TH2D '" << gen_hist_name << "' not found in " << gen_file << "\n";
    return 4;
  }
  if (hGen->GetSumw2N() == 0) hGen->Sumw2();

  // --- Read fit trees into maps keyed by (xq2bin, z_pt2_phi_bin) ---
  auto dataMap = ReadFitTreeToMap(&fData, fit_tree_name);
  auto mcMap   = ReadFitTreeToMap(&fMC,   fit_tree_name);

  std::cout << "Loaded DATA bins: " << dataMap.size() << "\n";
  std::cout << "Loaded MC   bins: " << mcMap.size()   << "\n";

  // --- Prepare output ---
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) { std::cerr << "ERROR: cannot create " << out_file << "\n"; return 5; }

  // Corrected: DATA * GEN / MC
  auto* hCorr = (TH2D*)hGen->Clone("h2_data_gen_over_mc");
  hCorr->Reset("ICES");
  hCorr->SetTitle("DATA * GEN / MC(rec);bin_xQ2;z_pt2_phi_bin");

  // Factor: GEN/MC(rec)
  auto* hGenOverMC = (TH2D*)hGen->Clone("h2_gen_over_mc");
  hGenOverMC->Reset("ICES");
  hGenOverMC->SetTitle("GEN / MC(rec);bin_xQ2;z_pt2_phi_bin");

  // NEW: Acceptance-like: MC(rec)/GEN (inverse of GEN/MC)
  auto* hMCOverGen = (TH2D*)hGen->Clone("h2_mc_over_gen");
  hMCOverGen->Reset("ICES");
  hMCOverGen->SetTitle("MC(rec) / GEN;bin_xQ2;z_pt2_phi_bin");

  // Optional but handy: write MC(rec) itself on the same 2D grid
  auto* hMCRec = (TH2D*)hGen->Clone("h2_mc_rec_input");
  hMCRec->Reset("ICES");
  hMCRec->SetTitle("MC(rec) from fits;bin_xQ2;z_pt2_phi_bin");

  // Output tree
  TTree tout("resp", "Per-bin DATA, MC(rec), GEN and combinations");

  int    out_xq2 = 0, out_z = 0;
  double out_data = 0, out_edata = 0;
  double out_mc   = 0, out_emc   = 0;
  double out_gen  = 0, out_egen  = 0;

  double out_corr = 0, out_ecorr = 0;
  double out_gen_over_mc = 0, out_egen_over_mc = 0;
  double out_mc_over_gen = 0, out_emc_over_gen = 0;   // NEW

  int out_status = 0;

  // status bits:
  //  1: missing DATA
  //  2: missing MC
  //  4: gen<=0
  //  8: mc<=0 (prevents division)
  tout.Branch("xq2bin", &out_xq2);
  tout.Branch("z_pt2_phi_bin", &out_z);

  tout.Branch("data", &out_data);
  tout.Branch("edata", &out_edata);

  tout.Branch("mc", &out_mc);
  tout.Branch("emc", &out_emc);

  tout.Branch("gen", &out_gen);
  tout.Branch("egen", &out_egen);

  tout.Branch("gen_over_mc", &out_gen_over_mc);
  tout.Branch("egen_over_mc", &out_egen_over_mc);

  tout.Branch("mc_over_gen", &out_mc_over_gen);       // NEW
  tout.Branch("emc_over_gen", &out_emc_over_gen);     // NEW

  tout.Branch("corr_data_gen_over_mc", &out_corr);
  tout.Branch("ecorr_data_gen_over_mc", &out_ecorr);

  tout.Branch("status", &out_status);

  // --- Loop over GEN histogram bins (defines the binning) ---
  const int nx = hGen->GetNbinsX();
  const int ny = hGen->GetNbinsY();

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {

      out_status = 0;

      // integerized bin centers (must match your encoding)
      out_xq2 = int(std::lround(hGen->GetXaxis()->GetBinCenter(ix)));
      out_z   = int(std::lround(hGen->GetYaxis()->GetBinCenter(iy)));

      // GEN from histogram
      out_gen  = hGen->GetBinContent(ix, iy);
      out_egen = hGen->GetBinError(ix, iy);
      if (out_gen <= 0) out_status |= 4;

      // DATA from map
      {
        auto it = dataMap.find(Key(out_xq2, out_z));
        if (it == dataMap.end()) {
          out_status |= 1;
          out_data = 0; out_edata = 0;
        } else {
          out_data  = it->second.v;
          out_edata = std::sqrt(std::max(0.0, it->second.e2));
        }
      }

      // MC(rec) from map
      {
        auto it = mcMap.find(Key(out_xq2, out_z));
        if (it == mcMap.end()) {
          out_status |= 2;
          out_mc = 0; out_emc = 0;
        } else {
          out_mc  = it->second.v;
          out_emc = std::sqrt(std::max(0.0, it->second.e2));
        }
      }
      if (out_mc <= 0) out_status |= 8;

      // Defaults
      out_gen_over_mc   = 0; out_egen_over_mc = 0;
      out_mc_over_gen   = 0; out_emc_over_gen = 0;
      out_corr          = 0; out_ecorr        = 0;

      // GEN/MC(rec) and DATA * GEN / MC(rec)
      if (out_gen > 0 && out_mc > 0) {
        out_gen_over_mc = out_gen / out_mc;

        // error for GEN / MC
        double rel2_f = 0.0;
        rel2_f += (out_egen / out_gen) * (out_egen / out_gen);
        rel2_f += (out_emc  / out_mc ) * (out_emc  / out_mc );
        out_egen_over_mc = std::abs(out_gen_over_mc) * std::sqrt(rel2_f);

        out_corr = out_data * out_gen_over_mc;

        // error for DATA * GEN / MC
        double rel2 = rel2_f;
        if (out_data > 0) rel2 += (out_edata / out_data) * (out_edata / out_data);
        out_ecorr = std::abs(out_corr) * std::sqrt(rel2);
      }

      // NEW: MC(rec) / GEN (acceptance-like)
      if (out_gen > 0) {
        out_mc_over_gen = out_mc / out_gen;

        // error for MC / GEN (only if both positive)
        if (out_mc > 0) {
          double rel2_a = 0.0;
          rel2_a += (out_emc  / out_mc ) * (out_emc  / out_mc );
          rel2_a += (out_egen / out_gen) * (out_egen / out_gen);
          out_emc_over_gen = std::abs(out_mc_over_gen) * std::sqrt(rel2_a);
        } else {
          out_emc_over_gen = 0.0;
        }
      }

      // Fill output histograms on the GEN grid
      hGenOverMC->SetBinContent(ix, iy, out_gen_over_mc);
      hGenOverMC->SetBinError(ix, iy, out_egen_over_mc);

      hMCOverGen->SetBinContent(ix, iy, out_mc_over_gen);
      hMCOverGen->SetBinError(ix, iy, out_emc_over_gen);

      hCorr->SetBinContent(ix, iy, out_corr);
      hCorr->SetBinError(ix, iy, out_ecorr);

      hMCRec->SetBinContent(ix, iy, out_mc);
      hMCRec->SetBinError(ix, iy, out_emc);

      tout.Fill();
    }
  }

  fout.cd();
  hGen->Write("h2_gen_input");
  hMCRec->Write();          // optional but useful
  hGenOverMC->Write();
  hMCOverGen->Write();      // NEW: rec/gen
  hCorr->Write();
  tout.Write();
  fout.Close();

  std::cout << "Wrote: " << out_file << "\n";
  return 0;
}
