// split_all_four_to_csv.cxx
// ROOT-loadable macro:
//   root -l
//   .L split_all_four_to_csv.cxx+
//   split_all_four_to_csv("measData_times_truth_over_meas.root","tables",8,11,8,true,-1);
//
// Writes 4 CSVs:
//   <prefix>_h_meas_mc.csv
//   <prefix>_h_true_mc.csv
//   <prefix>_h_meas_data.csv
//   <prefix>_h_measData_times_truth_over_meas.csv
//
// Columns: xq2_bin,ybin,z_bin,pt2_bin,phi_bin,phi_center_deg,content,error

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <memory>

#include "TFile.h"
#include "TH2.h"
#include "TAxis.h"

static inline int comp_index_1based(int z1b, int pt21b, int phi1b, int nPt2, int nPhi) {
  const int zpt2 = (z1b - 1) * nPt2 + pt21b;
  return (zpt2 - 1) * nPhi + phi1b; // 1..(nZ*nPt2*nPhi)
}

static inline double phi_center_deg(int phi1b, int nPhi) {
  const double w = 360.0 / (double)nPhi;
  return (phi1b - 0.5) * w;
}

static int infer_nPhi_from_h2(const TH2* h2, int nZ, int nPt2) {
  const int nx = h2->GetNbinsX();
  const int denom = nZ * nPt2;
  if (denom <= 0) return 0;

  // Common case from your producer: nx = denom*nPhi + 1
  if ((nx - 1) > 0 && ((nx - 1) % denom == 0)) return (nx - 1) / denom;

  // Alternate: nx = denom*nPhi
  if (nx % denom == 0) return nx / denom;

  return 0;
}

static bool export_th2_to_csv(const TH2* h2,
                              const std::string& csvPath,
                              int nZ, int nPt2, int nPhi,
                              bool skipZeros,
                              int maxXq2)
{
  std::ofstream out(csvPath);
  if (!out) {
    std::cerr << "ERROR: cannot open CSV for writing: " << csvPath << "\n";
    return false;
  }

  out << "xq2_bin,ybin,z_bin,pt2_bin,phi_bin,phi_center_deg,content,error\n";

  const TAxis* xax = h2->GetXaxis();
  const TAxis* yax = h2->GetYaxis();
  const int nx = h2->GetNbinsX();
  const int ny = h2->GetNbinsY();

  for (int iy = 1; iy <= ny; ++iy) {
    const double xq2_center = yax->GetBinCenter(iy);
    const int xq2_int = (int)std::llround(xq2_center);
    if (maxXq2 >= 0 && xq2_int > maxXq2) continue;

    for (int iz = 1; iz <= nZ; ++iz) {
      for (int ipt = 1; ipt <= nPt2; ++ipt) {
        for (int iphi = 1; iphi <= nPhi; ++iphi) {

          const int comp = comp_index_1based(iz, ipt, iphi, nPt2, nPhi);
          const int bx = xax->FindBin((double)comp);
          if (bx < 1 || bx > nx) continue;

          const double v = h2->GetBinContent(bx, iy);
          const double e = h2->GetBinError(bx, iy);

          if (skipZeros && v == 0.0 && e == 0.0) continue;

          out << xq2_int << ","
              << iy << ","
              << iz << ","
              << ipt << ","
              << iphi << ","
              << phi_center_deg(iphi, nPhi) << ","
              << v << ","
              << e << "\n";
        }
      }
    }
  }

  out.close();
  return true;
}

// ----------------------------------------------------------------------
// Entry point you call from ROOT
// ----------------------------------------------------------------------
void split_all_four_to_csv(const char* infile = "../measData_times_truth_over_meas.root",
                           const char* outPrefix = "tables",
                           int nZ = 8,
                           int nPt2 = 11,
                           int nPhi = 0,          // 0 => infer
                           bool skipZeros = false,
                           int maxXq2 = -1)       // -1 => all
{
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

  std::vector<TH2*> h2s;
  h2s.reserve(hnames.size());
  for (const auto& n : hnames) {
    auto* h = dynamic_cast<TH2*>(fin->Get(n.c_str()));
    if (!h) {
      std::cerr << "ERROR: TH2 '" << n << "' not found in " << infile << "\n";
      return;
    }
    h2s.push_back(h);
  }

  if (nPhi <= 0) {
    nPhi = infer_nPhi_from_h2(h2s[0], nZ, nPt2);
    if (nPhi <= 0) {
      std::cerr
        << "ERROR: could not infer nPhi from TH2 X bins.\n"
        << "Provide nPhi explicitly.\n"
        << "TH2 nbinsX=" << h2s[0]->GetNbinsX()
        << ", nZ*nPt2=" << (nZ * nPt2) << "\n";
      return;
    }
    std::cout << "Inferred nPhi = " << nPhi << "\n";
  }

  for (size_t i = 0; i < h2s.size(); ++i) {
    const std::string& hname = hnames[i];
    const TH2* h2 = h2s[i];

    const std::string csvPath = std::string(outPrefix) + "_" + hname + ".csv";
    std::cout << "Writing " << csvPath << " ...\n";

    if (!export_th2_to_csv(h2, csvPath, nZ, nPt2, nPhi, skipZeros, maxXq2)) {
      std::cerr << "ERROR exporting " << hname << "\n";
      return;
    }
  }

  fin->Close();
  std::cout << "Done.\n";
}
