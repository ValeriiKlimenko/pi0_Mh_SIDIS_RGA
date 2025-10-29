// File: define_bin_migr_hists_data.cxx
// Usage from ROOT:
//   root -l -q 'define_bin_migr_hists_data.cxx("path_to_folder_IN","unfolding","h22",80,80,120)'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH3.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <memory>


/*

[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005138_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005250_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005300_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005301_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005302_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005325_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005370_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005402_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005414_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005415_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005416_out.root
[SKIP] Missing required column in: /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/rec_data//data_f2018_005417_out.root

*/

// -------- Helpers

// Return a human-readable type for a branch: class name for objects, leaf type for PODs.
static std::string BranchTypeName(TTree* t, const std::string& col) {
  if (!t) return "";
  TBranch* br = t->GetBranch(col.c_str());        // non-const in ROOT 6.30
  if (!br) return "";
  const char* cls = br->GetClassName();
  if (cls && *cls) return std::string(cls);       // object branch (e.g. std::vector<float>)
  // POD / leaflist
  if (auto* leaves = br->GetListOfLeaves()) {
    if (leaves->GetSize() > 0) {
      if (auto* leaf = dynamic_cast<TLeaf*>(leaves->At(0))) {
        const char* ty = leaf->GetTypeName();
        if (ty) return std::string(ty);
      }
    }
  }
  return "";
}

static bool HasAllColumns(TTree* t, const std::vector<std::string>& cols) {
  if (!t) return false;
  for (const auto& c : cols) {
    if (!t->GetBranch(c.c_str())) return false;   // GetBranch is non-const
  }
  return true;
}

// List .root files in `dir` that contain TTree `tree` and where the columns exist
// with consistent types across all selected files.
static std::vector<std::string>
CollectFilesWithTreeAndColumns(const std::string& dir,
                               const std::string& tree,
                               const std::vector<std::string>& cols)
{
  std::vector<std::string> out;
  std::map<std::string,std::string> refType; // column -> type string from first good file

  TSystemDirectory sysDir(dir.c_str(), dir.c_str());
  TList* list = sysDir.GetListOfFiles();
  if (!list) {
    std::cerr << "[WARN] Directory not found or empty: " << dir << "\n";
    return out;
  }

  list->Sort();
  TIter next(list);
  while (TSystemFile* f = static_cast<TSystemFile*>(next())) {
    const char* nameC = f->GetName();
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

    TTree* tt = nullptr;
    tf->GetObject(tree.c_str(), tt);
    if (!tt) {
      std::cerr << "[SKIP] Tree '" << tree << "' not found in: " << fullpath << "\n";
      continue;
    }

    if (!HasAllColumns(tt, cols)) {
      std::cerr << "[SKIP] Missing required column in: " << fullpath << "\n";
      continue;
    }

    bool typesOK = true;
    for (const auto& c : cols) {
      std::string ty = BranchTypeName(tt, c);
      if (ty.empty()) { typesOK = false; std::cerr
        << "[SKIP] Could not determine type for '" << c << "' in " << fullpath << "\n"; break; }
      if (!refType.count(c)) refType[c] = ty;
      else if (refType[c] != ty) {
        typesOK = false;
        std::cerr << "[SKIP] Type mismatch for '" << c << "' in " << fullpath
                  << " : found '" << ty << "', expected '" << refType[c] << "'\n";
        break;
      }
    }
    if (typesOK) out.emplace_back(fullpath);
  }

  std::cout << "[INFO] Selected " << out.size()
            << " files with uniform types for tree '" << tree
            << "' from directory: " << dir << "\n";
  return out;
}

// -------- Main callable

void define_bin_migr_hists_data(const std::string& dir,
                                const std::string& outDir,
                                const std::string& tree    = "h22",
                                int nx                     = 80,
                                int ny                     = 80,
                                int nz                     = 40)
{
  // Unbuffered stdout so we see progress as it happens
  std::cout.setf(std::ios::unitbuf);

  //ROOT::EnableImplicitMT(3);

  // Axes/columns:
  const std::string aCol = "bin_xBQ2_Valerii";   // X-axis (integer bin id)
  const std::string xCol = "zpt2phit_8x8x9";     // Y-axis (integer/ID-like)
  const std::string zCol = "pi0_m";              // Z-axis (mass)

  // --- Build a schema-uniform file list
  const auto files = CollectFilesWithTreeAndColumns(dir, tree, {aCol, xCol, zCol});
  if (files.empty()) {
    std::cerr << "[ERROR] No input files with compatible schema were found in: " << dir << "\n";
    return;
  }

  // Input
  ROOT::RDataFrame df(tree, files);

  // -------- Ranges (lazy, single pass when RunGraphs is called)
  // Let RDF deduce types; don't force <int>.
  auto minA = df.Min(aCol);
  auto maxA = df.Max(aCol);
  auto minX = df.Min(xCol);
  auto maxX = df.Max(xCol);

  // Launch together, catch early if a file still misbehaves.
  try {
    ROOT::RDF::RunGraphs({minA, maxA, minX, maxX});
  } catch (const std::exception& e) {
    std::cerr << "[ERROR] While scanning min/max: " << e.what() << "\n";
    // Pinpoint the offending file
    TChain ch(tree.c_str());
    for (const auto& f : files) ch.Add(f.c_str());
    TTreeReader r(&ch);
    // Use generic reader values (choose the widest integer that fits many cases)
    // If your columns are floating point, change to double.
    // Safer approach: use TLeaf to check dynamically; here we just walk to find where it dies.
    for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
      auto st = r.SetEntry(i);
      if (st != TTreeReader::kEntryValid) {
        std::cerr << "[ERROR] Reader fails at global entry " << i
                  << " in file " << ch.GetFile()->GetName()
                  << " (status=" << int(st) << ")\n";
        break;
      }
    }
    return;
  }

  // Treat A/X as integer bin IDs -> +/- 0.5 around min/max
  double aMin = static_cast<double>(*minA) - 0.5;
  double aMax = static_cast<double>(*maxA) + 0.5;
  double xMin = static_cast<double>(*minX) - 0.5;
  double xMax = static_cast<double>(*maxX) + 0.5;

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
  // nz stays as passed-in (default 40 or your call)

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
