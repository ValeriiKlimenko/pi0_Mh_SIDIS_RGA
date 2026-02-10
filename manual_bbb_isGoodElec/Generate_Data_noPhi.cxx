// File: make_meas_data_only.cxx
//
// DATA-only version of your make_scaled_hist.cxx:
//   - reads the fitted DATA tree (h22_fit)
//   - fills ONLY h_meas_data (content = nPions, error = sqrt(sum(errPions^2)))
//   - supports "no-phi" combined binning: fewer X-axis bins on zpt2phit
//
// Usage from ROOT:
//   root -l -q 'make_meas_data_only.cxx("INPUT_fitted.root","OUT.root",1)'
//     where isNoPhi=1 => uses RESP nZ = N_Zbins * N_pTbins_with_overflow * 1 + 1
//           isNoPhi=0 => uses RESP nZ = N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1
//
// Notes:
//   - The branch name in the fitted DATA tree is still "z_pt2_phi_bin"
//     (even when isNoPhi==1, it just contains the no-phi combined IDs).
//   - Y axis is xq2bin with 21 bins (0..20) like before.

#include <iostream>
#include <memory>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TError.h"
#include "TH2D.h"

#include <ROOT/RDataFrame.hxx>

// Your binning
#include "binning_params.cxx"   // N_Zbins, N_pTbins_with_overflow, N_phiTrbins

static int make_meas_data_only(const char* in_fitted_file,
                              const char* out_file = "meas_data_only.root",
                              int isNoPhi = 1,
                              const char* treename = "h22_fit")
{
  if (!in_fitted_file || !*in_fitted_file) {
    std::cerr << "ERROR: input file path is empty\n";
    return 1;
  }
  if (gSystem->AccessPathName(in_fitted_file)) {
    std::cerr << "ERROR: input file not found: " << in_fitted_file << "\n";
    return 2;
  }

  // Check schema
  {
    TFile f(in_fitted_file, "READ");
    if (f.IsZombie()) {
      std::cerr << "ERROR: cannot open input: " << in_fitted_file << "\n";
      return 3;
    }
    TTree* t = nullptr;
    f.GetObject(treename, t);
    if (!t) {
      std::cerr << "ERROR: tree '" << treename << "' not found in " << in_fitted_file << "\n";
      return 4;
    }
    const bool ok =
      t->GetBranch("xq2bin") &&
      t->GetBranch("z_pt2_phi_bin") &&
      t->GetBranch("nPions") &&
      t->GetBranch("errPions");
    if (!ok) {
      std::cerr << "ERROR: missing one or more required branches in " << in_fitted_file << "\n"
                << "  need: xq2bin, z_pt2_phi_bin, nPions, errPions\n";
      return 5;
    }
  }

  // Disable MT because we fill ROOT histograms inside Foreach (not thread-safe)
  ROOT::DisableImplicitMT();

  // ---------------- Fixed geometry ----------------
  const int    nX   = 21;     // xq2bin: 0..20
  const double x_lo = -0.5;
  const double x_hi = x_lo + nX;

  const int nphi_eff = (isNoPhi == 1) ? 1 : N_phiTrbins;

  // Combined-bin axis size (fewer bins if no-phi)
  int    nZ   = N_Zbins * N_pTbins_with_overflow * nphi_eff + 1;
  double z_lo = -0.5;
  double z_hi = z_lo + nZ;

  std::cout << "[INFO] isNoPhi=" << isNoPhi
            << "  nphi_eff=" << nphi_eff
            << "  nZ=" << nZ << " (range " << z_lo << " .. " << z_hi << ")\n";

  // ---------------- Book histograms ----------------
  auto h_meas_data = std::make_unique<TH2D>(
      "h_meas_data",
      (isNoPhi ? "Measured data;zpt2phit_8x8x9_nophi;xq2bin"
               : "Measured data;zpt2phit_8x8x9;xq2bin"),
      nZ, z_lo, z_hi,
      nX, x_lo, x_hi
  );
  h_meas_data->Sumw2(false);          // we do custom errors from errPions
  h_meas_data->SetDirectory(nullptr);

  // accumulator for Σ(errPions^2) per bin
  auto h_err2 = std::make_unique<TH2D>(
      "h_meas_data_err2_acc",
      "err^2 accumulator;zpt2phit;xq2bin",
      nZ, z_lo, z_hi,
      nX, x_lo, x_hi
  );
  h_err2->SetDirectory(nullptr);

  // ---------------- Fill from RDataFrame ----------------
  ROOT::RDataFrame df(treename, in_fitted_file);

  auto good = df.Filter(
      [](double w, double ew) {
        return std::isfinite(w) && w > 0.0 && std::isfinite(ew) && ew >= 0.0;
      },
      {"nPions", "errPions"}
  );

  auto in_range = [](double v, double lo, double hi) {
    return std::isfinite(v) && v >= lo && v < hi;
  };

  long long skipped = 0;

  good.Foreach(
      [&](int zbin, int xbin, double wgt, double egt) {
        const double z = static_cast<double>(zbin);
        const double x = static_cast<double>(xbin);

        if (in_range(z, z_lo, z_hi) && in_range(x, x_lo, x_hi)) {
          h_meas_data->Fill(z, x, wgt);
          h_err2->Fill(z, x, egt * egt);
        } else {
          ++skipped;
        }
      },
      {"z_pt2_phi_bin", "xq2bin", "nPions", "errPions"}
  );

  if (skipped > 0) {
    Warning("make_meas_data_only",
            "Skipped %lld out-of-range entries (check bin-id ranges vs RESP geometry).",
            (long long)skipped);
  }

  // finalize bin errors: err = sqrt(sum(errPions^2))
  for (int ix = 1; ix <= h_meas_data->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h_meas_data->GetNbinsY(); ++iy) {
      const double e2 = h_err2->GetBinContent(ix, iy);
      h_meas_data->SetBinError(ix, iy, (e2 > 0.0 ? std::sqrt(e2) : 0.0));
    }
  }

  std::cout << "[INFO] Filled h_meas_data. Entries=" << h_meas_data->GetEntries()
            << "  Integral=" << h_meas_data->Integral() << "\n";

  // ---------------- Write output ----------------
  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: cannot create output file " << out_file << "\n";
    return 6;
  }

  h_meas_data->Write("h_meas_data");
  // optional: keep err2 acc for debugging (comment out if you don’t want it)
  h_err2->Write("h_meas_data_err2_acc");

  fout.Write();
  fout.Close();

  std::cout << "Wrote: " << out_file << " (hist: h_meas_data)\n";
  return 0;
}

// ROOT macro entry-point convenience (so you can call without naming the function)
int make_meas_data_only_driver(const char* in_fitted_file,
                              const char* out_file = "meas_data_only.root",
                              int isNoPhi = 1)
{
  return make_meas_data_only(in_fitted_file, out_file, isNoPhi, "h22_fit");
}
