// File: z_pt2_compare_from_dir.C
// Usage:
//   root -l -b -q 'z_pt2_compare_from_dir.C()'
//   root -l -b -q 'z_pt2_compare_from_dir.C("/path/to/dir","h22_z","out.root","out.pdf")'

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

namespace fs = std::filesystem;

// -------------------- binning constants (as in your code) --------------------
namespace bins {
  constexpr int N_xq2bins = 16 + 3;                 // 19 xQ2 bins
  constexpr int N_Zbins = 8;
  constexpr int N_pTbins = 10;
  constexpr int N_pTbins_with_overflow = N_pTbins + 1; // 11
  constexpr int N_phiTrbins = 9;

  inline const std::vector<double>& z_edges() {
    static const std::vector<double> z = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0};
    return z;
  }
  inline const std::vector<double>& pt2_edges() {
    // 11 bins (includes the "overflow" step as an explicit last edge)
    static const std::vector<double> p = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1.0,1.5};
    return p;
  }
}

// -------------------- helpers --------------------
inline double center_from_edges(const std::vector<double>& e, int bin1based) {
  if (bin1based < 1 || bin1based >= (int)e.size()) return -999.0;
  const int i0 = bin1based - 1;
  return 0.5*(e[i0] + e[i0+1]);
}

// Decode composite index (1-based) -> (z_bin, pt2_bin, phi_bin)
// Layout matches your earlier code: fastest index = phi, then pt2(with overflow), then z
inline bool decode_composite_1based(int comp, int& z_bin, int& pt2_bin, int& phi_bin) {
  using namespace bins;
  if (comp < 1) return false;

  phi_bin = ((comp - 1) % N_phiTrbins) + 1;

  const int zpt2_bin = ((comp - 1) / N_phiTrbins) + 1;
  const int per_z = N_pTbins_with_overflow; // 11

  z_bin   = ((zpt2_bin - 1) / per_z) + 1;   // 1..8
  pt2_bin = ((zpt2_bin - 1) % per_z) + 1;   // 1..11

  if (z_bin < 1 || z_bin > N_Zbins) return false;
  if (pt2_bin < 1 || pt2_bin > per_z) return false;
  if (phi_bin < 1 || phi_bin > N_phiTrbins) return false;
  return true;
}

// -------------------- main --------------------
void check_z_pt2(const char* dir_path = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/testing_binning/",
                            const char* tree_name = "h22_z",
                            const char* out_root  = "zpt2_compare.root",
                            const char* out_pdf   = "zpt2_compare.pdf")
{
  using namespace bins;

  gStyle->SetOptStat(0);

  // Prepare histograms per xQ2 bin
  const auto& Z = z_edges();
  const auto& P = pt2_edges();
  const int NZ = (int)Z.size() - 1;
  const int NP = (int)P.size() - 1;

  // ROOT wants non-const pointers for variable binning constructors
  Double_t* zEdges = const_cast<Double_t*>(Z.data());
  Double_t* pEdges = const_cast<Double_t*>(P.data());

  std::vector<TH2D*> h_true(N_xq2bins + 1, nullptr); // index 1..N_xq2bins
  std::vector<TH2D*> h_reco(N_xq2bins + 1, nullptr);

  for (int b = 1; b <= N_xq2bins; ++b) {
    h_true[b] = new TH2D(Form("h_true_xq2_%02d", b),
                         Form("TRUE: z vs #it{p}_{T}^{2}  (xQ2 bin %d);z;#it{p}_{T}^{2} [GeV^{2}]", b),
                         NZ, zEdges, NP, pEdges);
    h_reco[b] = new TH2D(Form("h_reco_xq2_%02d", b),
                         Form("RECO: z(center) vs #it{p}_{T}^{2}(center)  (xQ2 bin %d);z_{reco};#it{p}_{T,reco}^{2} [GeV^{2}]", b),
                         NZ, zEdges, NP, pEdges);
  }

  // Iterate all ROOT files in directory
  size_t files_ok = 0, files_bad = 0;
  for (const auto& ent : fs::directory_iterator(dir_path)) {
    if (!ent.is_regular_file()) continue;
    if (ent.path().extension() != ".root") continue;

    const std::string fpath = ent.path().string();
    TFile fin(fpath.c_str(), "READ");
    if (fin.IsZombie()) { std::cerr << "WARN: cannot open " << fpath << "\n"; ++files_bad; continue; }

    TTree* t = dynamic_cast<TTree*>(fin.Get(tree_name));
    if (!t) { std::cerr << "WARN: tree '" << tree_name << "' not found in " << fpath << "\n"; ++files_bad; continue; }

    auto has = [&](const char* br){ return t->GetBranch(br) != nullptr; };

    // Required branches for TRUE plot
    if (!has("bin_xBQ2_Valerii") || !has("z") || !has("pi0_sidis_PT2")) {
      std::cerr << "WARN: missing required branches in " << fpath
                << " (need bin_xBQ2_Valerii, z, pi0_sidis_PT2)\n";
      ++files_bad; continue;
    }

    // Try to find composite bin for reconstruction (prefer zpt2phit_8x8x9, else z_pt2_phi_bin)
    const bool has_compA = has("zpt2phit_8x8x9");
    const bool has_compB = has("z_pt2_phi_bin");
    const bool can_reco  = has_compA || has_compB;

    Int_t    xq2bin = 0;
    Double_t z_true = 0.0, pt2_true = 0.0;
    Int_t    comp   = 0;

    t->SetBranchAddress("bin_xBQ2_Valerii", &xq2bin);
    t->SetBranchAddress("z",                &z_true);
    t->SetBranchAddress("pi0_sidis_PT2",    &pt2_true);
    if (can_reco) {
      t->SetBranchAddress(has_compA ? "zpt2phit_8x8x9" : "z_pt2_phi_bin", &comp);
    }

    const Long64_t nent = t->GetEntries();
    for (Long64_t ie = 0; ie < nent; ++ie) {
      t->GetEntry(ie);
      if (xq2bin < 1 || xq2bin > N_xq2bins) continue;

      // TRUE fill
      h_true[xq2bin]->Fill(z_true, pt2_true);

      // RECO fill (derived from composite index -> bin centers)
      if (can_reco && comp > 0) {
        int z_bin = 0, pt2_bin = 0, phi_bin = 0;
        if (decode_composite_1based(comp, z_bin, pt2_bin, phi_bin)) {
          const double z_reco   = center_from_edges(Z, z_bin);
          const double pt2_reco = center_from_edges(P, pt2_bin);
          h_reco[xq2bin]->Fill(z_reco, pt2_reco);
        }
      }
    }

    fin.Close();
    ++files_ok;
  }

  // Save histograms and also make a side-by-side PDF
  TFile fout(out_root, "RECREATE");
  for (int b = 1; b <= N_xq2bins; ++b) {
    if (h_true[b]) h_true[b]->Write();
    if (h_reco[b]) h_reco[b]->Write();
  }
  fout.Close();

  // Multi-page PDF
  TCanvas c("c", "z vs pt2: TRUE vs RECO", 1200, 500);
  bool first = true;
  for (int b = 1; b <= N_xq2bins; ++b) {
    c.Clear();
    c.Divide(2,1);

    c.cd(1);
    if (h_true[b]) { h_true[b]->Draw("COLZ"); }

    c.cd(2);
    if (h_reco[b]) { h_reco[b]->Draw("COLZ"); }

    if (first) { c.Print((std::string(out_pdf) + "(").c_str()); first = false; }
    else       { c.Print(out_pdf); }
  }
  if (!first) c.Print((std::string(out_pdf) + ")").c_str()); // close PDF

  std::cout << "Processed files: ok=" << files_ok << ", skipped=" << files_bad << "\n"
            << "Wrote histograms to: " << out_root << "\n"
            << "Wrote plots (TRUE vs RECO per xQ2 bin) to: " << out_pdf << std::endl;
}
