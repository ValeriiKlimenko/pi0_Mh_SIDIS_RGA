// make_gen_binning_2D.cxx
// Build one 2D hist from all gen_f2018_*.root files under miss_dir.
// Skips broken/unusable ROOT files and reports how many were skipped.

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <ROOT/RDataFrame.hxx>

// If you already include this elsewhere in your build, keep it consistent.
#include "binning_params.cxx"   // defines N_Zbins, N_pTbins_with_overflow, N_phiTrbins

namespace RESP {
  static const int    nX   = 21;   // xq2bin: 0..20  (edges -0.5 .. 20.5)
  static const double x_lo = -0.5;
  static const double x_hi = x_lo + nX;

  static const double z_lo = -0.5;
}

static bool IsGoodInputFile(const std::string& path,
                            const char* treename,
                            const char* xbranch,
                            const char* ybranch,
                            std::string* reason_out = nullptr)
{
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie()) {
    if (reason_out) *reason_out = "cannot open / zombie";
    return false;
  }

  // Your broken example: "has no keys"
  if (f->GetNkeys() <= 0) {
    if (reason_out) *reason_out = "no keys in file";
    return false;
  }

  auto* t = dynamic_cast<TTree*>(f->Get(treename));
  if (!t) {
    if (reason_out) *reason_out = std::string("missing or non-TTree object: ") + treename;
    return false;
  }

  if (!t->GetBranch(xbranch)) {
    if (reason_out) *reason_out = std::string("missing branch: ") + xbranch;
    return false;
  }
  if (!t->GetBranch(ybranch)) {
    if (reason_out) *reason_out = std::string("missing branch: ") + ybranch;
    return false;
  }

  return true;
}

// Returns 0 on success
int make_gen_binning_2D(const char* miss_dir = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning",
                        const char* out_file = "gen_binning_2D.root",
                        int isNoPhi = 0)
{
  const char* treename = "h22";
  const char* xvar     = "bin_xBQ2_Valerii";
  const char* zvar     = (isNoPhi == 1) ? "zpt2phit_8x8x9_nophi" : "zpt2phit_8x8x9";

  // (1) Build + FILTER the file list single-threaded
  std::vector<std::string> files;
  int n_broken = 0;
  std::vector<std::string> broken_list;

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

      const std::string path = std::string(miss_dir) + "/" + nm.Data();

      std::string reason;
      if (!IsGoodInputFile(path, treename, xvar, zvar, &reason)) {
        ++n_broken;
        broken_list.push_back(path + "  (" + reason + ")");
        continue; // <-- skip broken
      }

      files.emplace_back(path);
    }
  }

  std::cout << "Good files: " << files.size()
            << "   |   Skipped broken/unusable: " << n_broken << "\n";

  if (n_broken > 0) {
    std::cerr << "Example broken files (up to 10 shown):\n";
    for (size_t i = 0; i < broken_list.size() && i < 10; ++i) {
      std::cerr << "  " << broken_list[i] << "\n";
    }
  }

  if (files.empty()) {
    std::cerr << "No usable input files in: " << miss_dir
              << " (pattern: gen_f2018_*.root)\n";
    return 2;
  }

  // (2) Enable implicit MT before creating the RDataFrame
  ROOT::EnableImplicitMT();

  using namespace RESP;

  const int nZ = (isNoPhi == 1)
      ? (N_Zbins * N_pTbins_with_overflow + 1)
      : (N_Zbins * N_pTbins_with_overflow * N_phiTrbins + 1);

  const double z_lo = RESP::z_lo;
  const double z_hi = z_lo + nZ;

  std::string htitle = std::string("Truth-like;") + xvar + ";" + zvar;

  try {
    ROOT::RDataFrame df(treename, files);

    auto h2 = df.Histo2D(
        ROOT::RDF::TH2DModel(
            "h2_binX_vs_z",
            htitle.c_str(),
            RESP::nX, RESP::x_lo, RESP::x_hi,
            nZ, z_lo, z_hi),
        xvar, zvar);

    TFile fout(out_file, "RECREATE");
    if (fout.IsZombie()) {
      std::cerr << "ERROR: cannot create output file: " << out_file << "\n";
      return 3;
    }

    h2->Write("h2_binX_vs_z");   // triggers event loop
    fout.Write();
    fout.Close();

    std::cout << "Wrote: " << out_file << " with histogram 'h2_binX_vs_z'\n";
    std::cout << "Skipped broken/unusable files total: " << n_broken << "\n";
    return 0;
  }
  catch (const std::exception& e) {
    std::cerr << "ERROR: RDataFrame failed: " << e.what() << "\n";
    std::cerr << "Skipped broken/unusable files total: " << n_broken << "\n";
    return 4;
  }
}