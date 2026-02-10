// split_all_four_to_csv_nophi.cxx
// ROOT-loadable macro (NO-PHI version):
//   root -l
//   .L split_all_four_to_csv_nophi.cxx+
//   split_all_four_to_csv_nophi("meas_data_only_nophi.root","tables_nophi",8,11,true,-1);
//
// Writes CSVs (only those that exist in the ROOT file):
//   <prefix>_h_meas_mc.csv
//   <prefix>_h_true_mc.csv
//   <prefix>_h_meas_data.csv
//   <prefix>_h_measData_times_truth_over_meas.csv
//
// Columns (no phi):
//   xq2_bin,ybin,z_bin,pt2_bin,content,error
//
// IMPORTANT:
//  - In "no-phi" scheme, the combined bin id is: comp = (zpt2-1)*1 + 1  => comp == zpt2
//  - So TH2 X-axis is effectively zpt2 combined index (1..nZ*nPt2), often plus +1 extra bin.
//  - This code finds the X bin via axis FindBin(comp), and also tries FindBin(comp+1) if needed.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <memory>

#include "TFile.h"
#include "TH2.h"
#include "TAxis.h"

// z1b in [1..nZ], pt21b in [1..nPt2]
static inline int comp_index_nophi_1based(int z1b, int pt21b, int nPt2) {
  return (z1b - 1) * nPt2 + pt21b; // 1..(nZ*nPt2)
}

// infer whether TH2 has "+1" extra X-bin pattern:
//   nbinsX == nZ*nPt2 + 1  (common with your producers)
static bool has_plus1_xbin(const TH2* h2, int nZ, int nPt2) {
  if (!h2) return false;
  const int nx = h2->GetNbinsX();
  return (nx == (nZ * nPt2 + 1));
}

// find x-bin robustly for comp in no-phi case
static int find_xbin_for_comp(const TH2* h2, int comp, int nZ, int nPt2) {
  const TAxis* xax = h2->GetXaxis();
  const int nx = h2->GetNbinsX();
  if (!xax) return -1;

  // First try direct
  int bx = xax->FindBin((double)comp);
  if (bx >= 1 && bx <= nx) return bx;

  // If histogram has "+1" convention, sometimes content is shifted by 1.
  // Try comp+1 as fallback.
  if (has_plus1_xbin(h2, nZ, nPt2)) {
    bx = xax->FindBin((double)(comp + 1));
    if (bx >= 1 && bx <= nx) return bx;
  }

  return -1;
}

static bool export_th2_to_csv_nophi(const TH2* h2,
                                   const std::string& csvPath,
                                   int nZ, int nPt2,
                                   bool skipZeros,
                                   int maxXq2)
{
  std::ofstream out(csvPath);
  if (!out) {
    std::cerr << "ERROR: cannot open CSV for writing: " << csvPath << "\n";
    return false;
  }

  out << "xq2_bin,ybin,z_bin,pt2_bin,content,error\n";

  const TAxis* xax = h2->GetXaxis();
  const TAxis* yax = h2->GetYaxis();
  if (!xax || !yax) {
    std::cerr << "ERROR: TH2 missing axes\n";
    return false;
  }

  const int ny = h2->GetNbinsY();

  for (int iy = 1; iy <= ny; ++iy) {
    const double xq2_center = yax->GetBinCenter(iy);
    const int xq2_int = (int)std::llround(xq2_center);
    if (maxXq2 >= 0 && xq2_int > maxXq2) continue;

    for (int iz = 1; iz <= nZ; ++iz) {
      for (int ipt = 1; ipt <= nPt2; ++ipt) {

        const int comp = comp_index_nophi_1based(iz, ipt, nPt2);
        const int bx = find_xbin_for_comp(h2, comp, nZ, nPt2);
        if (bx < 1) continue;

        const double v = h2->GetBinContent(bx, iy);
        const double e = h2->GetBinError(bx, iy);

        if (skipZeros && v == 0.0 && e == 0.0) continue;

        out << xq2_int << ","
            << iy << ","
            << iz << ","
            << ipt << ","
            << v << ","
            << e << "\n";
      }
    }
  }

  out.close();
  return true;
}

// ----------------------------------------------------------------------
// Entry point you call from ROOT (NO-PHI)
// ----------------------------------------------------------------------
void split_all_four_to_csv_nophi(const char* infile = "../meas_data_only.root",
                                const char* outPrefix = "tables_nophi",
                                int nZ = 8,
                                int nPt2 = 11,
                                bool skipZeros = false,
                                int maxXq2 = -1) // -1 => all
{
  // Try these in order; only export those found
  const std::vector<std::string> hnames = {
    "h_meas_mc",
    "h_true_mc",
    "h_meas_data",
    "h_measData_times_truth_over_meas"
  };

  std::unique_ptr<TFile> fin(TFile::Open(infile, "READ"));
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << infile << "\n";
    return;
  }

  int exported = 0;

  for (const auto& hname : hnames) {
    auto* h = dynamic_cast<TH2*>(fin->Get(hname.c_str()));
    if (!h) {
      std::cout << "[INFO] TH2 '" << hname << "' not found in " << infile << " (skipping)\n";
      continue;
    }

    const std::string csvPath = std::string(outPrefix) + "_" + hname + ".csv";
    std::cout << "Writing " << csvPath << " ...\n";

    if (!export_th2_to_csv_nophi(h, csvPath, nZ, nPt2, skipZeros, maxXq2)) {
      std::cerr << "ERROR exporting " << hname << "\n";
      return;
    }

    ++exported;
  }

  fin->Close();
  std::cout << "Done. Exported " << exported << " CSV file(s).\n";
}
