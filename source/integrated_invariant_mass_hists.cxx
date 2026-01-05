// g++ -O2 -std=c++17 `root-config --cflags --libs` make_hists_selected_dirs.cxx -o make_hists_sel
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <filesystem>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

namespace fs = std::filesystem;

static bool ends_with(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

static std::vector<std::string> collect_out_files_in_named_dirs(const fs::path& base_dir,
                                                                const std::string& target_dirname) {
  std::vector<std::string> files;
  if (!fs::exists(base_dir)) return files;

  for (auto const& entry : fs::recursive_directory_iterator(base_dir)) {
    if (!entry.is_regular_file()) continue;

    const fs::path p = entry.path();
    if (p.extension() != ".root") continue;
    if (!ends_with(p.filename().string(), "_out.root")) continue;

    const std::string parent = p.parent_path().filename().string();
    if (parent != target_dirname) continue;

    files.push_back(p.string());
  }

  std::sort(files.begin(), files.end());
  return files;
}

static bool tree_has_branch(const std::string& file, const char* tree, const char* branch) {
  TFile f(file.c_str(), "READ");
  if (f.IsZombie()) return false;
  auto* t = dynamic_cast<TTree*>(f.Get(tree));
  if (!t) return false;
  return t->GetListOfBranches() && t->GetListOfBranches()->FindObject(branch);
}

static void make_hists_for_dir(const fs::path& base_dir,
                               const std::string& dirname,
                               const fs::path& out_dir) {
  const char* tree_name = "h22";

  auto files = collect_out_files_in_named_dirs(base_dir, dirname);
  if (files.empty()) {
    std::cout << "[SKIP] No '*_out.root' found under '" << base_dir
              << "' in directories named '" << dirname << "'\n";
    return;
  }

  const std::string first = files.front();

  // Need pi0_m only for reconstructed/data directories
  const bool need_pi0m = (dirname == "rec_data" || dirname == "rec_true");

  // Branch checks
  if (!tree_has_branch(first, tree_name, "bin_xBQ2_Valerii") ||
      !tree_has_branch(first, tree_name, "zpt2phit_8x8x9") ||
      (need_pi0m && !tree_has_branch(first, tree_name, "pi0_m"))) {
    std::cerr << "[ERROR] Missing required branches in: " << first
              << " (dir=" << dirname << ")\n";
    return;
  }

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df(tree_name, files);

  const auto n = df.Count().GetValue();
  if (n == 0) {
    std::cout << "[SKIP] Zero entries for '" << dirname << "'\n";
    return;
  }

  // integer-like bin columns -> build binning from min/max
  const int x_min = (int)df.Min("bin_xBQ2_Valerii").GetValue();
  const int x_max = (int)df.Max("bin_xBQ2_Valerii").GetValue();
  const int y_min = (int)df.Min("zpt2phit_8x8x9").GetValue();
  const int y_max = (int)df.Max("zpt2phit_8x8x9").GetValue();

  const int n_x = std::max(1, x_max - x_min + 1);
  const int n_y = std::max(1, y_max - y_min + 1);

  // pi0_m binning (only used for rec/data 3D + optional 1D)
  const int    n_m   = 200;
  const double m_min = 0.0;
  const double m_max = 0.3;

  fs::create_directories(out_dir);
  const fs::path out_path = out_dir / (dirname + "_hists.root");

  std::cout << "[MAKE] " << dirname << ": " << files.size()
            << " files, " << n << " entries -> " << out_path << "\n";

  TFile fout(out_path.string().c_str(), "RECREATE");

  // ---------- 1D (keep as before) ----------
  auto h_x = df.Histo1D(
      {("h_" + dirname + "_bin_xBQ2_Valerii").c_str(),
       "bin_xBQ2_Valerii;bin_xBQ2_Valerii;Counts",
       n_x, x_min - 0.5, x_max + 0.5},
      "bin_xBQ2_Valerii");

  auto h_y = df.Histo1D(
      {("h_" + dirname + "_zpt2phit_8x8x9").c_str(),
       "zpt2phit_8x8x9;zpt2phit_8x8x9;Counts",
       n_y, y_min - 0.5, y_max + 0.5},
      "zpt2phit_8x8x9");

  h_x->Write();
  h_y->Write();

  // ---------- NEW: 2D or 3D ----------
  if (need_pi0m) {
    // 3D: X=bin_xBQ2_Valerii, Y=zpt2phit_8x8x9, Z=pi0_m
    auto h3 = df.Histo3D(
        {("h3_" + dirname + "_xq2_zpt2phi_pi0m").c_str(),
         "pi0_m vs xQ2-bin vs zpt2phi;bin_xBQ2_Valerii;zpt2phit_8x8x9;pi0_m (GeV)",
         n_x, x_min - 0.5, x_max + 0.5,
         n_y, y_min - 0.5, y_max + 0.5,
         n_m, m_min, m_max},
        "bin_xBQ2_Valerii", "zpt2phit_8x8x9", "pi0_m");

    h3->Write();

    // (optional, kept) 1D pi0_m
    auto h_m = df.Histo1D(
        {("h_" + dirname + "_pi0_m").c_str(),
         "pi0_m;pi0_m (GeV);Counts",
         n_m, m_min, m_max},
        "pi0_m");
    h_m->Write();

  } else {
    // gen_binning: 2D only (no pi0_m)
    auto h2 = df.Histo2D(
        {("h2_" + dirname + "_xq2_zpt2phi").c_str(),
         "counts vs xQ2-bin vs zpt2phi;bin_xBQ2_Valerii;zpt2phit_8x8x9;Counts",
         n_x, x_min - 0.5, x_max + 0.5,
         n_y, y_min - 0.5, y_max + 0.5},
        "bin_xBQ2_Valerii", "zpt2phit_8x8x9");

    h2->Write();
  }

  fout.Write();
  fout.Close();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr
      << "Usage:\n"
      << "  ./make_hists_sel <base_dir> [out_dir]\n\n"
      << "Creates one ROOT file per directory name: rec_data, rec_true, gen_binning\n";
    return 2;
  }

  const fs::path base_dir = argv[1];
  const fs::path out_dir  = (argc >= 3) ? fs::path(argv[2]) : (base_dir / "hists_selected");

  make_hists_for_dir(base_dir, "rec_data",    out_dir);
  make_hists_for_dir(base_dir, "rec_true",    out_dir);
  make_hists_for_dir(base_dir, "gen_binning", out_dir);

  return 0;
}

int run_make_hists_sel(const char* base_dir = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0",
                       const char* out_dir = nullptr)
{
  if (!base_dir || std::string(base_dir).empty()) {
    std::cerr << "[ERROR] run_make_hists_sel: base_dir is empty\n";
    return 2;
  }

  fs::path base = fs::path(base_dir);
  fs::path out  = (out_dir && std::string(out_dir).size())
                    ? fs::path(out_dir)
                    : (base / "hists_selected");

  make_hists_for_dir(base, "rec_data",    out);
  make_hists_for_dir(base, "rec_true",    out);
  make_hists_for_dir(base, "gen_binning", out);

  return 0;
}
