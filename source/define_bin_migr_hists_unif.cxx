// define_bin_migr_hists_unified.cxx
//
// Unified wrapper around the two "define_bin_migr_hists" producers.
// - Simulation mode  : builds TH3 per (A,B) and writes per-A files (original behavior).
// - Data mode        : builds a single TH3 (A, X, Z) and writes one file (original behavior).
//
// Entry points:
//   void define_bin_migr_hists( ... )          // simulation (kept for compatibility)
//   void define_bin_migr_hists_data( ... )     // data       (kept for compatibility)
//   void define_bin_migr_hists_switch(dir, outDir, tree="h22", nx=80, ny=80, nz=60,
//                                      const std::string& logic="sim"|"data")
//
// Example:
//   root -l -b -q \
//     -e '.L source/define_bin_migr_hists_unified.cxx+' \
//     -e 'define_bin_migr_hists_switch("/path/in","/path/out","h22",80,80,60,"sim")'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>

#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TDirectory.h>

#include <algorithm>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <utility>
#include <memory>

// ---------------------------------------------
// Helpers (shared)
// ---------------------------------------------
static void EnableMT_default(int nthreads = 3) {
  // Unbuffer stdout so progress prints even if something dies mid-run.
  std::cout.setf(std::ios::unitbuf);
  ROOT::EnableImplicitMT(nthreads);
}

// ----- Data helper: list files that actually have the requested TTree -----
static std::vector<std::string> CollectFilesWithTree(const std::string &dir,
                                                     const std::string &tree) {
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

// ---------------------------------------------
// SIMULATION implementation (original behavior)
// ---------------------------------------------
void define_bin_migr_hists(const std::string& dir,
                           const std::string& outDir,
                           const std::string& tree    = "h22",
                           int nx                     = 80,
                           int ny                     = 80,
                           int nz                     = 60) {
  EnableMT_default(3);

  // Column names
  const std::string xCol = "zpt2phit_8x8x9";
  const std::string yCol = "zpt2phit_8x8x9m";
  const std::string zCol = "pi0_m";
  const std::string aCol = "bin_xBQ2_Valerii";
  const std::string bCol = "bin_xBQ2_Valeriim";

  using AType = int;
  using BType = int;
  using Pair  = std::pair<AType,BType>;
  using PairSet = std::set<Pair>;

  ROOT::RDataFrame df(tree, dir + "/*.root");

  // -------- 1) Metadata: ranges + existing (A,B) pairs (lazy)
  auto minx  = df.Min<int>(xCol);
  auto maxx  = df.Max<int>(xCol);
  auto miny  = df.Min<int>(yCol);
  auto maxy  = df.Max<int>(yCol);

  auto dfPairs = df.Define("ab_pair",
                           [](AType a, BType b){ return Pair{a,b}; },
                           {aCol, bCol});

  // Aggregate into a set of pairs
  auto update_set = [](const PairSet &state, const Pair &p){
    PairSet out = state;
    out.insert(p);
    return out;
  };
  auto merge_sets = [](const PairSet &a, const PairSet &b){
    if (a.size() >= b.size()) { PairSet out = a; out.insert(b.begin(), b.end()); return out; }
    else                      { PairSet out = b; out.insert(a.begin(), a.end()); return out; }
  };

  const PairSet empty_set{};
  auto uniqPairsRes = dfPairs.Aggregate(update_set, merge_sets, "ab_pair", empty_set);

  // Run metadata actions in one pass
  ROOT::RDF::RunGraphs({minx, maxx, miny, maxy, uniqPairsRes});

  double xMin = *minx - 0.5, xMax = *maxx + 0.5;
  double yMin = *miny - 0.5, yMax = *maxy + 0.5;
  double zMin = 0.0,         zMax = 0.4;

  if (xMax <= xMin) xMax = xMin + 1.0;
  if (yMax <= yMin) yMax = yMin + 1.0;
  if (zMax <= zMin) zMax = zMin + 1.0;

  // Natural binning for integer-like axes
  nx = int(xMax - xMin);
  ny = int(yMax - yMin);

  PairSet uniqPairs = *uniqPairsRes;
  std::set<AType> uniqA;
  for (const auto &p : uniqPairs) uniqA.insert(p.first);

  std::cout << "Found " << uniqPairs.size() << " (A,B) pairs across "
            << uniqA.size() << " unique A values.\n";
  std::cout << "Binning: nx=" << nx << " ny=" << ny << " nz=" << nz
            << " | x:[" << xMin << "," << xMax << "]"
            << " y:[" << yMin << "," << yMax << "]"
            << " z:[" << zMin << "," << zMax << "]\n";

  // -------- 2) Book all histograms for existing pairs, single RunGraphs
  std::vector<ROOT::RDF::RResultHandle> actions;

  struct Pack {
    AType a;
    BType b;
    ROOT::RDF::RResultPtr<TH3D> h3;
  };
  std::map<AType, std::vector<Pack>> packsByA;

  auto book_hist_for = [&](AType aVal, BType bVal){
    auto dfAB = df.Filter([=](AType a, BType b){ return a==aVal && b==bVal; }, {aCol, bCol});

    std::ostringstream hname, htitle;
    hname  << "h3_" << aCol << "_" << aVal << "__" << bCol << "_" << bVal;
    htitle << "TH3 for " << aCol << "=" << aVal << ", " << bCol << "=" << bVal
           << "; " << xCol << "; " << yCol << "; " << zCol;

    ROOT::RDF::TH3DModel model(hname.str().c_str(), htitle.str().c_str(),
                               nx, xMin, xMax, ny, yMin, yMax, nz, zMin, zMax);

    auto h3 = dfAB.Histo3D(model, xCol, yCol, zCol);
    actions.emplace_back(h3);
    packsByA[aVal].push_back(Pack{aVal, bVal, h3});
  };

  for (const auto& pr : uniqPairs) book_hist_for(pr.first, pr.second);

  // Progress probe
  auto progress = df.Count();
  progress.OnPartialResult(1'000'000, [](ULong64_t &n){
    std::cout << "[RDF] processed " << n << " entries...\n";
  });
  actions.emplace_back(progress);

  std::cout << "Filling...\n";
  try {
    ROOT::RDF::RunGraphs(actions);
  } catch (const std::exception &e) {
    std::cerr << "RunGraphs threw: " << e.what() << std::endl;
    return;
  } catch (...) {
    std::cerr << "RunGraphs threw an unknown exception." << std::endl;
    return;
  }
  std::cout << "Done filling.\n";

  // --- Summary
  size_t nBooked = 0, nGT0 = 0, nGT10 = 0;
  for (auto &kv : packsByA) {
    nBooked += kv.second.size();
    for (auto &p : kv.second) {
      const auto entries = p.h3->GetEntries();
      if (entries > 0)  ++nGT0;
      if (entries > 10) ++nGT10;
    }
  }
  std::cout << "[Summary] Booked " << nBooked
            << " histograms; " << nGT0 << " have >0 entries; "
            << nGT10 << " have >10 entries.\n";

  // -------- 3) Write: only histograms with > 10 entries, grouped per A
  gSystem->MakeDirectory(outDir.c_str());

  size_t totalWritten = 0;
  for (auto &kv : packsByA) {
    const AType aVal = kv.first;
    std::cout << "Writing A=" << aVal << " ...\n";

    std::ostringstream fn;
    fn << outDir << "/h3_" << aCol << "_" << aVal << ".root";
    std::unique_ptr<TFile> fout(TFile::Open(fn.str().c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
      std::cerr << "ERROR: cannot create " << fn.str() << "\n";
      continue;
    }

    TDirectory* aDir = fout->mkdir((aCol + "_" + std::to_string(aVal)).c_str());
    aDir->cd();

    int written = 0;
    for (auto &p : kv.second) {
      if (p.h3->GetEntries() <= 10) continue; // keep only >10
      p.h3->Write();
      ++written;
      ++totalWritten;
      std::cout << "  Wrote " << p.h3->GetName()
                << " (entries=" << p.h3->GetEntries() << ")\n";
    }

    fout->Write();
    fout->Close();
    std::cout << "A=" << aVal << ": wrote " << written << " histograms (>10 entries).\n";
  }

  std::cout << "[Summary] Wrote " << totalWritten
            << " histograms with >10 entries in total.\n";
  std::cout << "Done. Output files are in: " << outDir << "\n";
}

// ---------------------------------------------
// DATA implementation (original behavior)
// ---------------------------------------------
void define_bin_migr_hists_data(const std::string& dir,
                                const std::string& outDir,
                                const std::string& tree    = "h22",
                                int nx                     = 80,
                                int ny                     = 80,
                                int nz                     = 60) {
  EnableMT_default(3);

  // Axes/columns:
  const std::string aCol = "bin_xBQ2_Valerii";   // X-axis (integer bin id)
  const std::string xCol = "zpt2phit_8x8x9";     // Y-axis (integer/ID-like)
  const std::string zCol = "pi0_m";              // Z-axis (mass)

  // Build the file list, skipping files missing the tree
  const auto files = CollectFilesWithTree(dir, tree);
  if (files.empty()) {
    std::cerr << "[ERROR] No input files with tree '" << tree << "' were found in: " << dir << "\n";
    return;
  }

  ROOT::RDataFrame df(tree, files);

  // -------- Ranges (lazy; computed on RunGraphs)
  auto minA = df.Min<int>(aCol);
  auto maxA = df.Max<int>(aCol);
  auto minX = df.Min<int>(xCol);
  auto maxX = df.Max<int>(xCol);

  ROOT::RDF::RunGraphs({minA, maxA, minX, maxX});

  // Treat A/X as integer bin IDs -> +/- 0.5 around min/max
  double aMin = *minA - 0.5, aMax = *maxA + 0.5;
  double xMin = *minX - 0.5, xMax = *maxX + 0.5;

  if (aMax <= aMin) aMax = aMin + 1.0;
  if (xMax <= xMin) xMax = xMin + 1.0;

  double zMin = 0.0, zMax = 0.4;
  if (zMax <= zMin) zMax = zMin + 1.0;

  // Natural binning for integer-like axes
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

  // Progress probe
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

  h3->Write();
  fout->Write();
  fout->Close();

  std::cout << "Wrote: " << fn.str() << "\n";
}

// ---------------------------------------------
// Unified switch entry
// ---------------------------------------------
void define_bin_migr_hists_switch(const std::string& dir,
                                  const std::string& outDir,
                                  const std::string& tree    = "h22",
                                  int nx                     = 80,
                                  int ny                     = 80,
                                  int nz                     = 40,
                                  const std::string& logic   = "sim") {
  if (logic == "data" || logic == "DATA") {
    define_bin_migr_hists_data(dir, outDir, tree, nx, ny, nz);
  } else if (logic == "sim" || logic == "SIM") {
    define_bin_migr_hists(dir, outDir, tree, nx, ny, nz);
  } else {
    std::cerr << "[ERROR] Unknown logic flag: " << logic << " (use \"sim\" or \"data\")\n";
  }
}
