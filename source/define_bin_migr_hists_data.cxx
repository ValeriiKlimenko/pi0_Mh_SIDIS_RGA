// File: define_bin_migr_hists.cxx
// Usage from ROOT:
//   root -l -q 'define_bin_migr_hists.cxx("path_to_folder_IN","unfolding","h22",80,80,120)'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TH3.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

// Helper: list .root files in `dir` that actually contain a TTree named `tree`
static std::vector<std::string> CollectFilesWithTree(const std::string &dir, const std::string &tree)
{
  std::vector<std::string> out;
  TSystemDirectory sysDir(dir.c_str(), dir.c_str());
  TList *list = sysDir.GetListOfFiles();
  if (!list) {
    std::cerr << "[WARN] Directory not found or empty: " << dir << "\n";
    return out;
  }

  list->Sort();
  TIter next(list);
  while (TSystemFile *f = static_cast<TSystemFile*>(next())) {
    const char *nameC = f->GetName();
    if (!nameC) continue;
    std::string name(nameC);
    if (f->IsDirectory()) continue;
    if (name == "." || name == "..") continue;
    if (name.size() < 6 || name.rfind(".root") != name.size()-5) continue;

    const std::string fullpath = dir + "/" + name;

    std::unique_ptr<TFile> tf(TFile::Open(fullpath.c_str(), "READ"));
    if (!tf || tf->IsZombie()) {
      std::cerr << "[SKIP] Cannot open file: " << fullpath << "\n";
      continue;
    }

    TTree *tt = nullptr;
    tf->GetObject(tree.c_str(), tt);
    if (!tt) {
      std::cerr << "[SKIP] Tree '" << tree << "' not found in: " << fullpath << "\n";
      continue;
    }

    out.emplace_back(fullpath);
  }

  std::cout << "[INFO] Selected " << out.size() << " files with tree '" << tree
            << "' from directory: " << dir << "\n";
  return out;
}

void define_bin_migr_hists_data(const std::string& dir,
                           const std::string& outDir,
                           const std::string& tree    = "h22",
                           int nx                     = 80,
                           int ny                     = 80,
                           int nz                     = 60)
{
  // Unbuffered stdout so we see progress even if something crashes mid-run
  std::cout.setf(std::ios::unitbuf);

  ROOT::EnableImplicitMT(3);

  // Axes/columns:
  const std::string aCol = "bin_xBQ2_Valerii";   // X-axis (integer bin id)
  const std::string xCol = "zpt2phit_8x8x9";     // Y-axis (integer/ID-like)
  const std::string zCol = "pi0_m";              // Z-axis (mass)

  // --- Build the file list, skipping files missing the tree
  const auto files = CollectFilesWithTree(dir, tree);
  if (files.empty()) {
    std::cerr << "[ERROR] No input files with tree '" << tree << "' were found in: " << dir << "\n";
    return;
  }

  // Input
  ROOT::RDataFrame df(tree, files);

  // -------- Ranges (lazy, single pass when RunGraphs is called)
  auto minA = df.Min<int>(aCol);
  auto maxA = df.Max<int>(aCol);
  auto minX = df.Min<int>(xCol);
  auto maxX = df.Max<int>(xCol);

  // Launch them together
  ROOT::RDF::RunGraphs({minA, maxA, minX, maxX});

  // Treat A/X as integer bin IDs -> +/- 0.5 around min/max
  double aMin = *minA - 0.5;
  double aMax = *maxA + 0.5;
  double xMin = *minX - 0.5;
  double xMax = *maxX + 0.5;

  // Guard against empty/degenerate ranges
  if (aMax <= aMin) aMax = aMin + 1.0;
  if (xMax <= xMin) xMax = xMin + 1.0;

  // Default z-range
  double zMin = 0.0;
  double zMax = 0.4;
  if (zMax <= zMin) zMax = zMin + 1.0;

  // If axes are integer IDs, use natural bin counts from ranges
  nx = static_cast<int>(aMax - aMin);
  ny = static_cast<int>(xMax - xMin);
  // nz stays as passed-in (default 60)

  std::cout << "Binning:\n"
            << "  A (" << aCol << "): nx=" << nx << " [" << aMin << "," << aMax << "]\n"
            << "  X (" << xCol << "): ny=" << ny << " [" << xMin << "," << xMax << "]\n"
            << "  Z (" << zCol << "): nz=" << nz << " [" << zMin << "," << zMax << "]\n";

  // -------- Book one TH3D over (A, X, Z)
  ROOT::RDF::TH3DModel model(
      "h3_binX_zpt_pi0m",
      (std::string("TH3; ") + aCol + ";" + xCol + ";" + zCol).c_str(),
      nx, aMin, aMax,
      ny, xMin, xMax,
      nz, zMin, zMax
  );

  auto h3 = df.Histo3D(model, aCol, xCol, zCol);

  // Progress probe during event loop
  auto progress = df.Count();
  progress.OnPartialResult(1'000'000, [](ULong64_t &n){
    std::cout << "[RDF] processed " << n << " entries...\n";
  });

  std::cout << "Filling...\n";
  try {
    ROOT::RDF::RunGraphs({h3, progress});
  } catch (const std::exception &e) {
    std::cerr << "RunGraphs threw: " << e.what() << std::endl;
    return;
  } catch (...) {
    std::cerr << "RunGraphs threw an unknown exception." << std::endl;
    return;
  }
  std::cout << "Done filling. Entries: " << h3->GetEntries() << "\n";

  // -------- Write a single output file
  gSystem->MakeDirectory(outDir.c_str());
  std::ostringstream fn;
  fn << outDir << "/h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m.root";

  std::unique_ptr<TFile> fout(TFile::Open(fn.str().c_str(), "RECREATE"));
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create " << fn.str() << "\n";
    return;
  }

  h3->Write();      // write the histogram
  fout->Write();
  fout->Close();

  std::cout << "Wrote: " << fn.str() << "\n";
}
