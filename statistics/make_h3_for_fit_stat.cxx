// make_h3_from_rdf.cxx
// Build TH3D (bin_xQ2 × zpt2phi × pi0_m) from all ROOT files in two folders (Data & MC)
// Outputs: h3_data.root, h3_mc.root in the same directory as this script.

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TFile.h>
#include <TH3D.h>
#include <TROOT.h>

#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

namespace fs = std::filesystem;

// ------------------- EDIT THESE TWO PATHS -------------------
static const char* DATA_DIR = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data/";  // folder with Data .root files
static const char* MC_DIR   = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/rec_true/";    // folder with MC .root files
// -----------------------------------------------------------

struct DatasetSpec {
  std::string name;   // "Data" or "MC"
  fs::path    dir;    // input folder
  bool        isMC;   // choose column names
};

static std::vector<std::string> list_root_files(const fs::path& p) {
  std::vector<std::string> out;
  if (!fs::exists(p) || !fs::is_directory(p)) return out;
  for (auto& e : fs::directory_iterator(p)) {
    if (e.is_regular_file()) {
      auto ext = e.path().extension().string();
      if (ext == ".root" || ext == ".ROOT") out.push_back(e.path().string());
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

static fs::path script_dir() {
  fs::path here = fs::path(__FILE__).parent_path();
  if (!here.empty() && fs::exists(here)) return here;
  return fs::current_path();
}

static void build_h3_for_dataset(const DatasetSpec& ds) {
  auto files = list_root_files(ds.dir);
  if (files.empty()) {
    std::cout << "[WARN] No .root files in " << ds.dir << " — skipping " << ds.name << "\n";
    return;
  }

  const char* treename   = "h22";  // your Snapshot tree name
  const char* col_xq2    = ds.isMC ? "bin_xBQ2_Valerii" : "bin_xBQ2_Valerii";
  // If you meant "bin_xBQ2_Valerii**m**" for MC, use that:
  col_xq2 = ds.isMC ? "bin_xBQ2_Valerii" : "bin_xBQ2_Valerii"; // adjust if needed

  const char* col_zpt2ph = ds.isMC ? "zpt2phit_8x8x9"   : "zpt2phit_8x8x9";
  const char* col_mass   = "pi0_m";
  //const char* col_inbins = ds.isMC ? "isEventINbins"   : "isEventINbins";

  ROOT::RDataFrame df(treename, files);

  auto& dfv = df;
  

  // Find max indices to size the TH3D so that bin centers are integers
  auto maxXv = dfv.Max<int>(col_xq2);
  auto maxYv = dfv.Max<int>(col_zpt2ph);
  const int maxX = static_cast<int>(*maxXv);
  const int maxY = static_cast<int>(*maxYv);

  if (maxX <= 0 || maxY <= 0) {
    std::cout << "[WARN] Non-positive maxX/maxY (" << maxX << "," << maxY
              << ") for " << ds.name << " — skipping\n";
    return;
  }

  ROOT::RDF::TH3DModel model(
    ds.isMC ? "h3_binX_zpt_pi0m_mc" : "h3_binX_zpt_pi0m_data",
    (std::string(";bin_xQ2 (") + ds.name + ");zpt2phi (" + ds.name + ");#pi^{0} mass (GeV)").c_str(),
    maxX, 0.5, maxX + 0.5,
    maxY, 0.5, maxY + 0.5,
    60,   0.0, 0.4   // adjust mass binning if needed
  );

  auto h3 = dfv.Histo3D(model, col_xq2, col_zpt2ph, col_mass);

  fs::path outdir = script_dir();
  fs::path out    = outdir / (ds.isMC ? "h3_mc.root" : "h3_data.root");

  TFile fout(out.string().c_str(), "RECREATE");
  h3->Write();   // name: model.GetName()
  fout.Close();

  std::cout << "[OK] Wrote " << out << "  ("
            << h3->GetName() << ": "
            << maxX << "×" << maxY << "×" << h3->GetNbinsZ()
            << " bins)\n";
}

// ROOT entry point with no parameters (uses the hard-coded paths)
void make_h3_for_fit_stat() {
  ROOT::EnableImplicitMT();
  build_h3_for_dataset(DatasetSpec{"Data", DATA_DIR, false});
  build_h3_for_dataset(DatasetSpec{"MC",   MC_DIR,   true });
}
