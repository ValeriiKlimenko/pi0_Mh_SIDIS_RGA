// File: z_pt2_compare_from_dir.C
// Usage:
//   root -l -b -q 'z_pt2_compare_from_dir.C()'
//   root -l -b -q 'z_pt2_compare_from_dir.C("/path/to/dir","h22","out.root","out.pdf")'

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
#include "TLegend.h"
#include "TLatex.h"

namespace fs = std::filesystem;

// -------------------- binning constants --------------------
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
// Layout: fastest index = phi(9), then pt2(with overflow=11), then z(8)
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

// Make a 0/1 "mask" (coverage) from a 2D histogram (any content -> 1)
static TH2D* make_mask(const TH2D* src, const char* name, const char* title) {
  auto* m = (TH2D*)src->Clone(name);
  m->SetTitle(title);
  m->Reset("ICES"); // keep axes, zero contents
  const int nx = src->GetNbinsX();
  const int ny = src->GetNbinsY();
  for (int ix=1; ix<=nx; ++ix) {
    for (int iy=1; iy<=ny; ++iy) {
      m->SetBinContent(ix, iy, src->GetBinContent(ix,iy) > 0 ? 1.0 : 0.0);
    }
  }
  return m;
}

// Compute RECO/TRUE ratio where TRUE>0; 0 otherwise
static TH2D* make_ratio(const TH2D* reco, const TH2D* tru, const char* name, const char* title) {
  auto* r = (TH2D*)tru->Clone(name);
  r->SetTitle(title);
  r->Reset("ICES");
  const int nx = tru->GetNbinsX();
  const int ny = tru->GetNbinsY();
  for (int ix=1; ix<=nx; ++ix) {
    for (int iy=1; iy<=ny; ++iy) {
      const double t = tru->GetBinContent(ix,iy);
      const double q = reco->GetBinContent(ix,iy);
      r->SetBinContent(ix, iy, (t>0.0) ? (q/t) : 0.0);
    }
  }
  return r;
}

// Build categorical coverage map:
// 0 = none, 1 = TRUE-only, 2 = RECO-only, 3 = both
static TH2D* make_covcat(const TH2D* mask_true, const TH2D* mask_reco,
                         const char* name, const char* title) {
  auto* c = (TH2D*)mask_true->Clone(name);
  c->SetTitle(title);
  c->Reset("ICES");
  const int nx = c->GetNbinsX();
  const int ny = c->GetNbinsY();
  for (int ix=1; ix<=nx; ++ix) {
    for (int iy=1; iy<=ny; ++iy) {
      const int mt = (mask_true->GetBinContent(ix,iy) > 0) ? 1 : 0;
      const int mr = (mask_reco->GetBinContent(ix,iy) > 0) ? 1 : 0;
      const int cat = mt + 2*mr; // 0,1,2,3 with meanings above
      c->SetBinContent(ix,iy,cat);
    }
  }
  c->SetMinimum(0.0);
  c->SetMaximum(3.0);
  return c;
}

// -------------------- main --------------------
void z_pt2_compare_from_dir(const char* dir_path = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/testing_binning/",
                            const char* tree_name = "h22",
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
                         Form("TRUE: z vs #it{p}_{T}^{2}  (xQ^{2} bin %d);z;#it{p}_{T}^{2} [GeV^{2}]", b),
                         NZ, zEdges, NP, pEdges);
    h_reco[b] = new TH2D(Form("h_reco_xq2_%02d", b),
                         Form("RECO: z(center) vs #it{p}_{T}^{2}(center)  (xQ^{2} bin %d);z_{reco};#it{p}_{T,reco}^{2} [GeV^{2}]", b),
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

    // Composite bin for reconstruction (prefer zpt2phit_8x8x9, else z_pt2_phi_bin)
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

  // ---- Coverage products per xQ2 bin
  std::vector<TH2D*> h_true_mask(N_xq2bins + 1, nullptr);
  std::vector<TH2D*> h_reco_mask(N_xq2bins + 1, nullptr);
  std::vector<TH2D*> h_covcat   (N_xq2bins + 1, nullptr); // 0 none, 1 TRUE-only, 2 RECO-only, 3 both
  std::vector<TH2D*> h_ratio    (N_xq2bins + 1, nullptr); // RECO/TRUE where TRUE>0

  // Summary tree (one row per xQ2,zBin,pT2Bin)
  TFile fout(out_root, "RECREATE");
  int   T_xq2=0, T_zBin=0, T_pT2Bin=0;
  int   T_trueAny=0, T_recoAny=0, T_cat=0;
  float T_trueCnt=0.f, T_recoCnt=0.f, T_ratio=0.f;
  TTree covTree("coverage", "Per-(xQ2,z,pT2) coverage and counts");
  covTree.Branch("xq2",     &T_xq2,     "xq2/I");
  covTree.Branch("zBin",    &T_zBin,    "zBin/I");
  covTree.Branch("pT2Bin",  &T_pT2Bin,  "pT2Bin/I");
  covTree.Branch("trueAny", &T_trueAny, "trueAny/I");
  covTree.Branch("recoAny", &T_recoAny, "recoAny/I");
  covTree.Branch("cat",     &T_cat,     "cat/I"); // 0 none,1 TRUE-only,2 RECO-only,3 both
  covTree.Branch("trueCnt", &T_trueCnt, "trueCnt/F");
  covTree.Branch("recoCnt", &T_recoCnt, "recoCnt/F");
  covTree.Branch("ratio",   &T_ratio,   "ratio/F"); // RECO/TRUE if TRUE>0 else 0

  for (int b = 1; b <= N_xq2bins; ++b) {
    h_true_mask[b] = make_mask(h_true[b], Form("h_true_mask_xq2_%02d", b),
                               Form("TRUE coverage (any>0)  xQ^{2} bin %d;z;#it{p}_{T}^{2} [GeV^{2}]", b));
    h_reco_mask[b] = make_mask(h_reco[b], Form("h_reco_mask_xq2_%02d", b),
                               Form("RECO coverage (any>0)  xQ^{2} bin %d;z;#it{p}_{T}^{2} [GeV^{2}]", b));
    h_covcat[b]    = make_covcat(h_true_mask[b], h_reco_mask[b],
                                 Form("h_covcat_xq2_%02d", b),
                                 Form("Coverage category (0 none,1 TRUE,2 RECO,3 both)  xQ^{2} bin %d;z;#it{p}_{T}^{2} [GeV^{2}]", b));
    h_ratio[b]     = make_ratio(h_reco[b], h_true[b],
                                Form("h_ratio_xq2_%02d", b),
                                Form("RECO/TRUE per bin (TRUE>0)  xQ^{2} bin %d;z;#it{p}_{T}^{2} [GeV^{2}]", b));

    // Write histos
    h_true[b]->Write();
    h_reco[b]->Write();
    h_true_mask[b]->Write();
    h_reco_mask[b]->Write();
    h_covcat[b]->Write();
    h_ratio[b]->Write();

    // Fill summary tree
    const int nx = h_true[b]->GetNbinsX();
    const int ny = h_true[b]->GetNbinsY();
    for (int ix=1; ix<=nx; ++ix) {
      for (int iy=1; iy<=ny; ++iy) {
        const double tcnt = h_true[b]->GetBinContent(ix,iy);
        const double rcnt = h_reco[b]->GetBinContent(ix,iy);
        const int tAny = tcnt>0 ? 1:0;
        const int rAny = rcnt>0 ? 1:0;
        const int cat  = tAny + 2*rAny;
        T_xq2=b; T_zBin=ix; T_pT2Bin=iy;
        T_trueAny=tAny; T_recoAny=rAny; T_cat=cat;
        T_trueCnt=tcnt; T_recoCnt=rcnt;
        T_ratio = (tcnt>0) ? float(rcnt/tcnt) : 0.f;
        covTree.Fill();
      }
    }

    // Simple text summary per xQ2
    int nb_true=0, nb_reco=0, nb_both=0, nb_tonly=0, nb_ronly=0;
    for (int ix=1; ix<=h_covcat[b]->GetNbinsX(); ++ix) {
      for (int iy=1; iy<=h_covcat[b]->GetNbinsY(); ++iy) {
        const int cat = (int)h_covcat[b]->GetBinContent(ix,iy);
        nb_true += ((h_true_mask[b]->GetBinContent(ix,iy)>0)?1:0);
        nb_reco += ((h_reco_mask[b]->GetBinContent(ix,iy)>0)?1:0);
        if      (cat==3) ++nb_both;
        else if (cat==1) ++nb_tonly;
        else if (cat==2) ++nb_ronly;
      }
    }
    std::cout << Form("[xQ2 bin %2d] covered TRUE bins=%3d, RECO bins=%3d, BOTH=%3d, TRUE-only=%3d, RECO-only=%3d\n",
                      b, nb_true, nb_reco, nb_both, nb_tonly, nb_ronly);
  }

  covTree.Write();

  for (int b = 1; b <= N_xq2bins; ++b) {
    h_true[b]->Write();  h_true[b]->SetDirectory(nullptr);
    h_reco[b]->Write();  h_reco[b]->SetDirectory(nullptr);
    h_true_mask[b]->Write(); h_true_mask[b]->SetDirectory(nullptr);
    h_reco_mask[b]->Write(); h_reco_mask[b]->SetDirectory(nullptr);
    h_covcat[b]->Write();    h_covcat[b]->SetDirectory(nullptr);
    h_ratio[b]->Write();     h_ratio[b]->SetDirectory(nullptr);
  }

  
  fout.Close();

  // ---- Multi-page PDF (4 panels per xQ2): TRUE, RECO, Coverage category, Ratio
  TCanvas c("c", "z vs pT^{2} coverage", 1600, 900);
  bool first = true;
  for (int b = 1; b <= N_xq2bins; ++b) {
    c.Clear();
    c.Divide(2,2);

    c.cd(1);
    h_true[b]->Draw("COLZ");
    TLatex l1; l1.SetNDC(); l1.SetTextSize(0.04);
    l1.DrawLatex(0.12,0.93,Form("xQ^{2} bin %d  — TRUE counts", b));

    c.cd(2);
    h_reco[b]->Draw("COLZ");
    TLatex l2; l2.SetNDC(); l2.SetTextSize(0.04);
    l2.DrawLatex(0.12,0.93,Form("xQ^{2} bin %d  — RECO counts", b));

    c.cd(3);
    h_covcat[b]->SetMinimum(0); h_covcat[b]->SetMaximum(3);
    h_covcat[b]->Draw("COLZ");
    TLatex lt; lt.SetNDC(); lt.SetTextSize(0.035);
    lt.DrawLatex(0.12,0.93,"Coverage category: 0 none, 1 TRUE-only, 2 RECO-only, 3 both");

    c.cd(4);
    h_ratio[b]->SetMinimum(0); h_ratio[b]->Draw("COLZ");
    TLatex lr; lr.SetNDC(); lr.SetTextSize(0.04);
    lr.DrawLatex(0.12,0.93,"RECO/TRUE (only where TRUE>0)");

    if (first) { c.Print((std::string(out_pdf) + "(").c_str()); first = false; }
    else       { c.Print(out_pdf); }
  }
  if (!first) c.Print((std::string(out_pdf) + ")").c_str()); // close PDF

  std::cout << "Processed files: ok=" << files_ok << ", skipped=" << files_bad << "\n"
            << "Wrote histograms and coverage maps to: " << out_root << "\n"
            << "Wrote plots to: " << out_pdf << std::endl;
}

// Optional: keep your old entry point name for convenience
void check_z_pt2(const char* dir_path = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/testing_binning/",
                 const char* tree_name = "h22",
                 const char* out_root  = "zpt2_compare.root",
                 const char* out_pdf   = "zpt2_compare.pdf") {
  z_pt2_compare_from_dir(dir_path, tree_name, out_root, out_pdf);
}
