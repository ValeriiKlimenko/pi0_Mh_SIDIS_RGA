// make_gen_binning_2D.cxx
// Build one 2D hist from all gen_f2018_*.root files under miss_dir.

#include <iostream>
#include <vector>
#include <string>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TFile.h>
#include <TH2D.h>
#include <ROOT/RDataFrame.hxx>


// If you already include this elsewhere in your build, keep it consistent.
#include "binning_params.cxx"   // defines N_Zbins, N_pTbins_with_overflow, N_phiTrbins

namespace RESP {
  static const int    nX   = 21;   // xq2bin: 0..20  (edges -0.5 .. 20.5)
  static const double x_lo = -0.5;
  static const double x_hi = x_lo + nX;

  static const int    nZ   = N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1;
  static const double z_lo = -0.5;
  static const double z_hi = z_lo + nZ;
}

// Returns 0 on success
int make_gen_binning_2D(const char* miss_dir = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning",
                        const char* out_file = "gen_binning_2D.root")
{
  // (1) Build the file list single-threaded as before
  std::vector<std::string> files;
  {
    TSystemDirectory dir("miss_dir", miss_dir);
    TList* l = dir.GetListOfFiles();
    if (!l) { std::cerr << "No directory listing for: " << miss_dir << "\n"; return 1; }

    TIter it(l);
    while (TSystemFile* f = (TSystemFile*)it()) {
      if (f->IsDirectory()) continue;
      TString nm = f->GetName();
      if (!nm.BeginsWith("gen_f2018_")) continue;
      if (!nm.EndsWith(".root", TString::kIgnoreCase)) continue;
      files.emplace_back(std::string(miss_dir) + "/" + nm.Data());
    }
  }
  if (files.empty()) {
    std::cerr << "No matching input files in: " << miss_dir
              << " (pattern: gen_f2018_*.root)\n";
    return 2;
  }

  // (2) Enable implicit MT before creating the RDataFrame
  //     You can pass a number to cap threads, e.g. EnableImplicitMT(8).
  ROOT::EnableImplicitMT();

  using namespace RESP;
  ROOT::RDataFrame df("h22", files);

  auto h2 = df.Histo2D(
      ROOT::RDF::TH2DModel(
          "h2_binX_vs_z",
          "Truth-like;bin_xBQ2_Valerii;zpt2phit_8x8x9",
          RESP::nX, RESP::x_lo, RESP::x_hi,
          RESP::nZ, RESP::z_lo, RESP::z_hi),
      "bin_xBQ2_Valerii", "zpt2phit_8x8x9");

  TFile fout(out_file, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: cannot create output file: " << out_file << "\n";
    return 3;
  }
  h2->Write("h2_binX_vs_z");   // triggers the event loop; thread-safe merge
  fout.Write();
  fout.Close();

  // Optional: turn MT back off for later single-threaded code
  // ROOT::DisableImplicitMT();

  std::cout << "Wrote: " << out_file << " with histogram 'h2_binX_vs_z'\n";
  return 0;
}
