// File: define_bin_migr_hists.cxx
// Usage from ROOT:
//   root -l -q 'define_bin_migr_hists.cxx("path_to_folder_IN","h22","unfolding",80,80,120)'

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
#include <TFile.h>
#include <TSystem.h>
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

// /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/rec_response/ no last two cuts (on gen lvl)

void define_bin_migr_hists(const std::string& dir,
                           const std::string& outDir,
                           const std::string& tree    = "h22",
                           int nx                     = 80,
                           int ny                     = 80,
                           int nz                     = 40)
{
  // Flush std::cout on each << so messages show even if something dies mid-run
  std::cout.setf(std::ios::unitbuf);

  //ROOT::EnableImplicitMT(3);

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

  // -------- 1) Metadata: ranges + existing (A,B) pairs (lazy & memory-light)
  auto minx  = df.Min<int>(xCol);
  auto maxx  = df.Max<int>(xCol);
  auto miny  = df.Min<int>(yCol);
  auto maxy  = df.Max<int>(yCol);

  // Single (A,B) column
  auto dfPairs = df.Define("ab_pair",
                           [](AType a, BType b){ return Pair{a,b}; },
                           {aCol, bCol});

  // ROOT 6.30 Aggregate requires: (state, value) -> new_state  and  (state, state) -> merged_state
  auto update_set = [](const PairSet &state, const Pair &p){
    PairSet out = state;
    out.insert(p);
    return out;
  };
  auto merge_sets = [](const PairSet &a, const PairSet &b){
    // insert smaller into larger
    if (a.size() >= b.size()) { PairSet out = a; out.insert(b.begin(), b.end()); return out; }
    else                      { PairSet out = b; out.insert(a.begin(), a.end()); return out; }
  };

  const PairSet empty_set{};
  auto uniqPairsRes = dfPairs.Aggregate(update_set, merge_sets, "ab_pair", empty_set);

  // Run metadata actions together (single pass)
  ROOT::RDF::RunGraphs({minx, maxx, miny, maxy, uniqPairsRes});

  double xMin = *minx - 0.5;
  double xMax = *maxx + 0.5;
  double yMin = *miny - 0.5;
  double yMax = *maxy + 0.5;
  double zMin = 0.0;
  double zMax = 0.4;

  if (xMax <= xMin) xMax = xMin + 1.0;
  if (yMax <= yMin) yMax = yMin + 1.0;
  if (zMax <= zMin) zMax = zMin + 1.0;

  // If x/y are integer bin IDs, keep auto-binning like before:
  nx = int(xMax - xMin);
  ny = int(yMax - yMin);

  PairSet uniqPairs = *uniqPairsRes;
  std::set<AType> uniqA;
  for (const auto &p : uniqPairs) uniqA.insert(p.first);

  std::cout << "Found " << uniqPairs.size() << " existing (A,B) pairs across "
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
    ROOT::RDF::RResultPtr<TH3D> h3;  // TH3DModel for ROOT 6.30 compatibility
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

  // Progress probe so you see output during the event loop
  auto progress = df.Count();
  progress.OnPartialResult(1'000'000, [](ULong64_t &n){
    std::cout << "[RDF] processed " << n << " entries...\n";
  });
  actions.emplace_back(progress);

  std::cout << "before filling\n";

  try {
    ROOT::RDF::RunGraphs(actions);
    std::cout << "after filling\n";
  } catch (const std::exception &e) {
    std::cerr << "RunGraphs threw: " << e.what() << std::endl;
    return;
  } catch (...) {
    std::cerr << "RunGraphs threw an unknown exception." << std::endl;
    return;
  }

  // --- Summary: how many histograms have entries?
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
      if (p.h3->GetEntries() <= 3) continue; // keep only >3
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
