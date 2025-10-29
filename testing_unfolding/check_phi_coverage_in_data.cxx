// File: add_zpt2phi_branches.C
// Usage:
//   root -l -b -q 'add_zpt2phi_branches.C("my_fitted.root")'
//   root -l -b -q 'add_zpt2phi_branches.C("my_fitted.root", "h22_fit", "my_fitted_with_angles.root")'

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

// -------------------- binning constants: MATCH zpt2phit_8x8x9 --------------------
namespace bins {
  constexpr int N_xq2bins = 16 + 3;   // context only
  constexpr int N_Zbins   = 8;        // 8 z bins
  constexpr int N_pTbins  = 11;       
  constexpr int N_phiTrbins = 8;      // 9 phi bins

  inline const std::vector<double>& z_edges() {
    // 8 bins => 9 edges
    static const std::vector<double> z = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0};
    return z;
  }
  inline const std::vector<double>& pt2_edges() {
    // Provide 8-bin pT^2 edges (adjust if your exact 8-bin scheme differs)
    // This choice merges the higher pT^2 region into one bin [0.5,0.8].
    static const std::vector<double> p = {0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.65,0.8,1,1.5};
    return p;
  }
}

// -------------------- helpers --------------------
inline double center_from_edges(const std::vector<double>& e, int bin1based) {
  if (bin1based < 1 || bin1based >= (int)e.size()) return -999.0;
  const int i0 = bin1based - 1;
  return 0.5*(e[i0] + e[i0+1]);
}

inline bool decode_composite_1based(int comp, int& z_bin, int& pt2_bin, int& phi_bin) {
  using namespace bins;
  if (comp < 1) return false;

  // Packing: phi is fastest; then pt2; then z
  phi_bin = ((comp - 1) % N_phiTrbins) + 1;
  const int zpt2_bin = ((comp - 1) / N_phiTrbins) + 1;
  const int per_z = N_pTbins; // 8

  z_bin   = ((zpt2_bin - 1) / per_z) + 1;   // 1..N_Zbins
  pt2_bin = ((zpt2_bin - 1) % per_z) + 1;   // 1..N_pTbins

  if (z_bin   < 1 || z_bin   > N_Zbins)     return false;
  if (pt2_bin < 1 || pt2_bin > per_z)       return false;
  if (phi_bin < 1 || phi_bin > N_phiTrbins) return false;
  return true;
}

inline double phi_center_deg(int phi_bin) {
  using namespace bins;
  const double w = 360.0 / N_phiTrbins;     // 40°
  return (phi_bin - 0.5)*w;                 // 20, 60, …, 340
}

// -------------------- main driver --------------------
void check_phi_coverage_in_data(const char* in_file = "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_rec_data/h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root",
                                const char* tree_name = "h22_fit",
                                const char* out_file  = "data_cover.root")
{
  using namespace bins;

  // Open input
  TFile fin(in_file, "READ");
  if (fin.IsZombie()) {
    throw std::runtime_error(std::string("Cannot open input ROOT file: ") + in_file);
  }
  TTree* tin = dynamic_cast<TTree*>(fin.Get(tree_name));
  if (!tin) {
    throw std::runtime_error(std::string("Tree not found: ") + tree_name);
  }

  // ----- bind the input branches (kept unchanged in output) -----
  Int_t    xq2bin = 0;
  Int_t    zpt2phi_comp = 0;        // "z_pt2_phi_bin"
  Int_t    zpt2phi_hist_bin = 0;    // "z_pt2_phi_hist_bin"
  Double_t nPions = 0, errPions = 0;

  tin->SetBranchAddress("xq2bin",             &xq2bin);
  tin->SetBranchAddress("z_pt2_phi_bin",      &zpt2phi_comp);
  tin->SetBranchAddress("z_pt2_phi_hist_bin", &zpt2phi_hist_bin);
  tin->SetBranchAddress("nPions",             &nPions);
  tin->SetBranchAddress("errPions",           &errPions);

  // Output file
  std::string outpath = (out_file && std::string(out_file).size())
                        ? std::string(out_file)
                        : (std::string(in_file).substr(0, std::string(in_file).find_last_of('.')) + "_with_angles.root");
  TFile fout(outpath.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output ROOT file: ") + outpath);
  }

  // Clone keeps ALL original branches
  TTree* tout = tin->CloneTree(0);

  // New branches: raw bin numbers + centers
  Int_t    z_bin = 0, pt2_bin = 0, phi_bin = 0;  // raw indices (1-based)
  Double_t z = 0.0, pt2 = 0.0, phi = 0.0;        // bin centers

  tout->Branch("z_bin",   &z_bin);
  tout->Branch("pt2_bin", &pt2_bin);
  tout->Branch("phi_bin", &phi_bin);

  tout->Branch("z",   &z);
  tout->Branch("pt2", &pt2);
  tout->Branch("phi", &phi);

  const auto& Z = z_edges();
  const auto& P = pt2_edges();

  // Sanity: edges must match declared binning
  if ((int)Z.size() != N_Zbins + 1) {
    throw std::runtime_error("z_edges size mismatch: expected " + std::to_string(N_Zbins+1));
  }
  if ((int)P.size() != N_pTbins + 1) {
    throw std::runtime_error("pt2_edges size mismatch: expected " + std::to_string(N_pTbins+1));
  }

  const int maxComp = N_phiTrbins * N_Zbins * N_pTbins; // 9*8*8 = 576

  // --- Counters
  Long64_t total_valid = 0, total_nonzero = 0;
  Long64_t total_failed = 0;
  std::vector<Long64_t> per_xq2_valid   (N_xq2bins + 1, 0);
  std::vector<Long64_t> per_xq2_nonzero (N_xq2bins + 1, 0);
  std::vector<Long64_t> per_xq2_failed  (N_xq2bins + 1, 0);
  Long64_t xq2_out_of_range_valid   = 0,
           xq2_out_of_range_nonzero = 0,
           xq2_out_of_range_failed  = 0;

  const Long64_t nent = tin->GetEntries();
  for (Long64_t ie = 0; ie < nent; ++ie) {
    tin->GetEntry(ie);

    // Prefer the raw Y-bin index if it’s in range; else use the composite value
    int comp = zpt2phi_comp;

    bool ok = decode_composite_1based(comp, z_bin, pt2_bin, phi_bin);
    if (ok) {
      z   = center_from_edges(Z, z_bin);
      pt2 = center_from_edges(P, pt2_bin);
      phi = phi_center_deg(phi_bin);

      ++total_valid;
      const bool nonzero = (nPions > 0.0);
      if (nonzero) ++total_nonzero;

      if (xq2bin >= 1 && xq2bin <= N_xq2bins) {
        ++per_xq2_valid[xq2bin];
        if (nonzero) ++per_xq2_nonzero[xq2bin];
      } else {
        ++xq2_out_of_range_valid;
        if (nonzero) ++xq2_out_of_range_nonzero;
      }
    } else {
      ++total_failed;
      if (xq2bin >= 1 && xq2bin <= N_xq2bins) {
        ++per_xq2_failed[xq2bin];
      } else {
        ++xq2_out_of_range_failed;
      }

      // mark invalid decode in output branches
      z_bin = pt2_bin = phi_bin = -1;
      z = pt2 = phi = -999.0;
    }

    // Fill copies all original branches (from CloneTree) + our new ones
    tout->Fill();
  }

  fout.cd();
  tout->Write();
  fout.Close();
  fin.Close();

  // -------------------- Summary printout --------------------
  std::cout << "\n=== Coverage summary for file: " << in_file << " ===\n";
  std::cout << "Entries read:                  " << nent << "\n";
  std::cout << "Valid (decodable) points:      " << total_valid << "\n";
  std::cout << "  of which nPions > 0:         " << total_nonzero << "\n";
  std::cout << "Failed to decode:              " << total_failed << "\n";
  std::cout << "Sanity (valid + failed):       " << (total_valid + total_failed) << "\n";

  std::cout << "\nPer xQ2 bin (1.." << N_xq2bins
            << ") => valid / (nPions>0) / failed\n";

  for (int b = 1; b <= N_xq2bins; ++b) {
    if (per_xq2_valid[b] == 0 && per_xq2_nonzero[b] == 0 && per_xq2_failed[b] == 0) continue;
    std::cout << "  xQ2 bin " << std::setw(2) << b << ": "
              << per_xq2_valid[b] << " / "
              << per_xq2_nonzero[b] << " / "
              << per_xq2_failed[b] << "\n";
  }
  if (xq2_out_of_range_valid > 0 || xq2_out_of_range_nonzero > 0 || xq2_out_of_range_failed > 0) {
    std::cout << "  xQ2 out-of-range: "
              << xq2_out_of_range_valid << " / "
              << xq2_out_of_range_nonzero << " / "
              << xq2_out_of_range_failed << "\n";
  }

  std::cout << "\nWrote tree with original branches + raw indices (z_bin, pt2_bin, phi_bin) "
               "and centers (z, pt2, phi) to:\n  "
            << outpath << "\n" << std::endl;
}
