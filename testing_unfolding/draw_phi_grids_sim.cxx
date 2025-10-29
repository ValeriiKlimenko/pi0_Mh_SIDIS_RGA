// draw_phi_grids_sim.cxx
// Make per-xQ2 canvases; grid of (z rows) x (pt2 columns); each cell: nPions(φ)
// Works on the SIMULATION outputs produced by split_and_fit_unified(..., Logic::Sim)
// Tree: h22_fit with branches: xq2bin, z_pt2_phi_bin, nPions, errPions (and *_gen variants)

#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <set>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TDirectory.h"

// -------------------- binning constants (SIM path; 8×11×8) --------------------
namespace bins {
  constexpr int N_Zbins      = 8;
  constexpr int N_pTbins     = 11;  // including your high-pT^2 tail bin
  constexpr int N_phi_bins   = 8;

  inline const std::vector<double>& z_edges() {
    static const std::vector<double> z = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0}; // 8 bins
    return z;
  }
  inline const std::vector<double>& pt2_edges() {
    // 11 bins (12 edges)
    static const std::vector<double> p = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1.0,1.5};
    return p;
  }

  inline double center_from_edges(const std::vector<double>& e, int bin1based) {
    if (bin1based < 1 || bin1based >= (int)e.size()) return -999.0;
    const int i0 = bin1based - 1;
    return 0.5*(e[i0] + e[i0+1]);
  }

  // composite decoding (1-based), packing: phi fastest, then pt2, then z
  inline bool decode_composite_1based(int comp, int& z_bin, int& pt2_bin, int& phi_bin) {
    if (comp < 1) return false;
    phi_bin = ((comp - 1) % N_phi_bins) + 1;
    const int zpt2_bin = ((comp - 1) / N_phi_bins) + 1;
    z_bin   = ((zpt2_bin - 1) / N_pTbins) + 1;
    pt2_bin = ((zpt2_bin - 1) % N_pTbins) + 1;
    if (z_bin   < 1 || z_bin   > N_Zbins)    return false;
    if (pt2_bin < 1 || pt2_bin > N_pTbins)   return false;
    if (phi_bin < 1 || phi_bin > N_phi_bins) return false;
    return true;
  }

  inline double phi_center_deg(int phi_bin) {
    const double w = 360.0 / N_phi_bins;       // 45°
    return (phi_bin - 0.5)*w;                  // 22.5, 67.5, …, 337.5
  }
}

// -------------------- aggregator --------------------
struct Cell {
  double y[bins::N_phi_bins]{};
  double ey2[bins::N_phi_bins]{};
  bool   has[bins::N_phi_bins]{};
};

using KeyXQ2      = int;
using KeyZPT2     = std::pair<int,int>;   // (z_bin, pt2_bin)
using GridForXQ2  = std::map<KeyZPT2, Cell>;
using AllData     = std::map<KeyXQ2, GridForXQ2>;

// -------------------- core --------------------
void draw_phi_grids_sim_impl(const std::vector<std::string>& files,
                             const std::string& outDir = "phi_grids",
                             bool useGen = false)
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir.c_str(), /*recursive=*/true);

  // ---- Read tree(s)
  TChain ch("h22_fit");
  for (const auto& f : files) ch.Add(f.c_str());

  int    xq2bin = 0, xq2bin_gen = 0;
  int    comp   = 0, comp_gen   = 0;
  double nPions = 0, errPions = 0;

  ch.SetBranchAddress("xq2bin", &xq2bin);
  if (ch.GetBranch("xq2bin_gen")) ch.SetBranchAddress("xq2bin_gen", &xq2bin_gen);
  ch.SetBranchAddress("z_pt2_phi_bin", &comp);
  if (ch.GetBranch("z_pt2_phi_bin_gen")) ch.SetBranchAddress("z_pt2_phi_bin_gen", &comp_gen);
  ch.SetBranchAddress("nPions", &nPions);
  ch.SetBranchAddress("errPions", &errPions);

  AllData agg;
  std::set<int> xq2_seen;

  const Long64_t nent = ch.GetEntries();
  for (Long64_t i=0;i<nent;i++) {
    ch.GetEntry(i);

    const int use_comp = useGen ? comp_gen : comp;
    int zbin=0, pt2bin=0, phibin=0;
    if (!bins::decode_composite_1based(use_comp, zbin, pt2bin, phibin)) continue;

    const int ixq2 = xq2bin;  // reconstructed xQ2 for counting/plotting
    xq2_seen.insert(ixq2);

    Cell& cell = agg[ixq2][{zbin, pt2bin}];
    const int iphi = phibin - 1;

    cell.y[iphi]   += nPions;            // sum yields
    cell.ey2[iphi] += errPions*errPions; // combine errors in quadrature
    cell.has[iphi]  = true;
  }

  // ---- For each xQ2 → make a grid canvas (rows=z, cols=pt2)
  const auto& zE  = bins::z_edges();
  const auto& pE  = bins::pt2_edges();

  for (int xq2 : xq2_seen) {
    // find y-max across all subpads for this xq2 to unify scales
    double yMax = 0.0;
    for (int z=1; z<=bins::N_Zbins; ++z) {
      for (int p=1; p<=bins::N_pTbins; ++p) {
        const auto it = agg[xq2].find({z,p});
        if (it == agg[xq2].end()) continue;
        const Cell& c = it->second;
        for (int ip=0; ip<bins::N_phi_bins; ++ip) {
          if (c.has[ip]) yMax = std::max(yMax, c.y[ip]);
        }
      }
    }
    if (yMax <= 0) yMax = 1.0;

    const int NX = bins::N_pTbins;
    const int NY = bins::N_Zbins;

    std::ostringstream ctitle;
    ctitle << "nPions vs #phi — xQ2 bin " << xq2
           << (useGen ? " (GEN z,pT^{2},#phi)" : " (RECO z,pT^{2},#phi)");

    TCanvas* c = new TCanvas(Form("c_phi_xq2_%d_%s", xq2, useGen?"gen":"reco"),
                             ctitle.str().c_str(), 2200, 1500);
    c->Divide(NX, NY, 0.0001, 0.0001);  // tight packing

    // Labels helpers
    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.05);

    for (int row=1; row<=NY; ++row) {
      const int zbin = row;          // top row = zbin 1 (you can flip if preferred)
      const double zc = bins::center_from_edges(zE, zbin);

      for (int col=1; col<=NX; ++col) {
        const int padIdx = (row-1)*NX + col;
        c->cd(padIdx);
        gPad->SetMargin(0.12, 0.02, 0.18, 0.05);

        const int pt2bin = col;
        const double pt2c = bins::center_from_edges(pE, pt2bin);

        // Frame (common axes)
        TH1F* frame = new TH1F(Form("frm_%d_%d_%d_%s", xq2, zbin, pt2bin, useGen?"g":"r"),
                               "", 100, 0.0, 360.0);
        frame->SetDirectory(nullptr);
        frame->SetMinimum(0.0);
        frame->SetMaximum(yMax*1.15);
        frame->GetXaxis()->SetTitle("#phi [deg]");
        frame->GetYaxis()->SetTitle("n_{#pi^{0}}");
        frame->GetXaxis()->SetLabelSize(0.06);
        frame->GetYaxis()->SetLabelSize(0.06);
        frame->GetXaxis()->SetTitleSize(0.07);
        frame->GetYaxis()->SetTitleSize(0.07);
        frame->GetYaxis()->SetTitleOffset(0.6);
        frame->GetXaxis()->SetTitleOffset(0.9);
        // show fewer axis titles to reduce clutter
        if (row != NY) frame->GetXaxis()->SetTitle("");
        if (col != 1 ) frame->GetYaxis()->SetTitle("");

        frame->Draw("AXIS");

        // Fetch cell data
        auto it = agg[xq2].find({zbin, pt2bin});
        if (it == agg[xq2].end()) {
          lat.SetTextSize(0.06);
          lat.DrawLatex(0.18, 0.80, Form("z=%.2f", zc));
          lat.DrawLatex(0.18, 0.68, Form("p_{T}^{2}=%.3g", pt2c));
          lat.SetTextSize(0.07);
          lat.DrawLatex(0.25, 0.45, "no data");
          continue;
        }

        const Cell& cell = it->second;

        // Build graph: only keep points that exist
        std::vector<double> vx, vy, vex, vey;
        vx.reserve(bins::N_phi_bins);
        vy.reserve(bins::N_phi_bins);
        vex.reserve(bins::N_phi_bins);
        vey.reserve(bins::N_phi_bins);

        for (int ip=0; ip<bins::N_phi_bins; ++ip) {
          if (!cell.has[ip]) continue;
          const int phi_bin = ip+1;
          vx.push_back(bins::phi_center_deg(phi_bin));
          vy.push_back(cell.y[ip]);
          vex.push_back(0.0);
          vey.push_back(std::sqrt(cell.ey2[ip]));
        }

        if (!vx.empty()) {
          TGraphErrors* gr = new TGraphErrors((int)vx.size(),
                                              vx.data(), vy.data(),
                                              vex.data(), vey.data());
          gr->SetMarkerStyle(20);
          gr->SetMarkerSize(0.9);
          gr->SetLineWidth(2);
          gr->Draw("P SAME");
        }

        // cell label (top-left)
        lat.SetTextSize(0.055);
        lat.DrawLatex(0.15, 0.86, Form("z=%.2f", zc));
        lat.DrawLatex(0.15, 0.74, Form("p_{T}^{2}=%.3g", pt2c));
      }
    }

    // Big title on the whole canvas
    c->cd(0);
    TLatex head; head.SetNDC(); head.SetTextAlign(13);
    head.SetTextSize(0.035);
    head.DrawLatex(0.02, 0.98, ctitle.str().c_str());

    // Save
    const std::string base = Form("%s/phi_grid_xq2_%d_%s", outDir.c_str(), xq2, useGen?"gen":"reco");
    c->SaveAs((base + ".png").c_str());
    c->SaveAs((base + ".pdf").c_str());
  }

  if (xq2_seen.empty()) {
    std::cerr << "[draw_phi_grids_sim] No entries found. Check input files/branches.\n";
  } else {
    std::cout << "[draw_phi_grids_sim] Wrote " << xq2_seen.size()
              << " canvases to: " << outDir << "\n";
  }
}

// -------------------- friendly entry point for ROOT -q --------------------
void draw_phi_grids_sim(const char* commaSeparatedFiles = "",
                        const char* outDir = "phi_grids",
                        bool useGen = false)
{
  // parse comma-separated list
  std::vector<std::string> files;
  {
    std::stringstream ss(commaSeparatedFiles ? commaSeparatedFiles : "");
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      if (!tok.empty()) files.push_back(tok);
    }
  }
  if (files.empty()) {
    std::cerr << "[draw_phi_grids_sim] Please pass at least one _fitted.root file.\n";
    return;
  }
  draw_phi_grids_sim_impl(files, outDir ? outDir : "phi_grids", useGen);
}
