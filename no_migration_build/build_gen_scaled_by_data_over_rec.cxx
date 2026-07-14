// build_gen_scaled_by_data_over_rec.cxx
//
// Creates output histogram/tree with, for each (ix,iy):
//   gen(ix,iy) * nPions_data(ix,iy) / nPions_rec(ix,iy)
//
// Inputs (hardcoded folder):
//   rec_data_hists_sliced_fitted.root  -> TTree "h22_fit"
//   rec_true_hists_sliced_fitted.root  -> TTree "h22_fit"   (this is your "rec")
//   gen_binning_hists.root             -> TH2D (bin_xBQ2_Valerii vs zpt2phit_8x8x9)
//
// Output:
//   gen_scaled_data_over_rec.root containing:
//     - TH2D h2_gen_scaled
//     - TH2D h2_scale_data_over_rec
//     - TTree scale_tree (one row per (ix,iy))
//
// Build:
//   g++ -O2 -std=c++17 `root-config --cflags --libs` build_gen_scaled_by_data_over_rec.cxx -o build_gen_scaled
//
// Run:
//   ./build_gen_scaled
//
// ROOT:
//   root -l
//   root [0] .L build_gen_scaled_by_data_over_rec.cxx+
//   root [1] run_build_gen_scaled();   // default paths

#include <map>
#include <utility>
#include <string>
#include <iostream>
#include <memory>
#include <cstring>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TClass.h"
#include "TH2D.h"

namespace {

// ---------- hardcoded paths ----------
static const char* kDir =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/";

static const char* kDataFitFile =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/rec_data_hists_sliced_fitted.root";

static const char* kRecFitFile =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/rec_true_hists_sliced_fitted.root";

static const char* kGenHistsFile =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/gen_binning_hists.root";

static const char* kDefaultOut =
  "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/hists_selected/gen_scaled_data_over_rec.root";

// ---------- helpers ----------
struct FitVal {
  double n   = -1.0;
  double err = -1.0;
  int    ok  = 0;
};

using Key = std::pair<int,int>; // (ix,iy)

static std::map<Key, FitVal> LoadFitMap(TTree* t, const char* tag) {
  std::map<Key, FitVal> m;
  if (!t) {
    std::cerr << "[ERROR] null tree for " << tag << "\n";
    return m;
  }

  int ix=0, iy=0, ok=0;
  double nPions=0, errPions=0;

  // required branches from your slicer:
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("ix", 1);
  t->SetBranchStatus("iy", 1);
  t->SetBranchStatus("ok", 1);
  t->SetBranchStatus("nPions", 1);
  t->SetBranchStatus("errPions", 1);

  t->SetBranchAddress("ix", &ix);
  t->SetBranchAddress("iy", &iy);
  t->SetBranchAddress("ok", &ok);
  t->SetBranchAddress("nPions", &nPions);
  t->SetBranchAddress("errPions", &errPions);

  const Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i) {
    t->GetEntry(i);

    // keep only good fits with positive yield
    if (ok != 1 || !(nPions > 0.0)) continue;

    Key k{ix,iy};
    FitVal v{nPions, errPions, ok};

    // if duplicates exist, keep the one with larger nPions (or smaller rel error)
    auto it = m.find(k);
    if (it == m.end()) {
      m.emplace(k, v);
    } else {
      const double rel_old = (it->second.n > 0) ? std::abs(it->second.err / it->second.n) : 1e9;
      const double rel_new = (v.n > 0) ? std::abs(v.err / v.n) : 1e9;
      if (v.n > it->second.n || rel_new < rel_old) it->second = v;
    }
  }

  std::cout << "[INFO] Loaded " << m.size() << " good fits from " << tag << "\n";
  return m;
}

static TH2D* FindFirstTH2D(TFile* f, const char* preferred_name = "h2_gen_binning_xq2_zpt2phi") {
  if (!f) return nullptr;

  if (preferred_name && std::strlen(preferred_name)) {
    if (auto* h = dynamic_cast<TH2D*>(f->Get(preferred_name))) return h;
  }

  // fallback: first TH2D at top-level
  TIter nextKey(f->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(nextKey())) {
    auto* cls = TClass::GetClass(key->GetClassName());
    if (!cls) continue;
    if (cls->InheritsFrom(TH2D::Class())) {
      return dynamic_cast<TH2D*>(key->ReadObj());
    }
  }
  return nullptr;
}

} // namespace

int BuildGenScaled(const char* out_path = nullptr) {
  const std::string outp = (out_path && std::strlen(out_path)) ? out_path : kDefaultOut;

  // open input files
  std::unique_ptr<TFile> fData(TFile::Open(kDataFitFile, "READ"));
  std::unique_ptr<TFile> fRec (TFile::Open(kRecFitFile,  "READ"));
  std::unique_ptr<TFile> fGen (TFile::Open(kGenHistsFile, "READ"));

  if (!fData || fData->IsZombie()) { std::cerr << "[ERROR] cannot open " << kDataFitFile << "\n"; return 2; }
  if (!fRec  || fRec->IsZombie())  { std::cerr << "[ERROR] cannot open " << kRecFitFile  << "\n"; return 2; }
  if (!fGen  || fGen->IsZombie())  { std::cerr << "[ERROR] cannot open " << kGenHistsFile << "\n"; return 2; }

  auto* tData = dynamic_cast<TTree*>(fData->Get("h22_fit"));
  auto* tRec  = dynamic_cast<TTree*>(fRec->Get("h22_fit"));
  if (!tData) { std::cerr << "[ERROR] missing TTree h22_fit in data file\n"; return 3; }
  if (!tRec)  { std::cerr << "[ERROR] missing TTree h22_fit in rec file\n";  return 3; }

  TH2D* hGen = FindFirstTH2D(fGen.get(), "h2_gen_binning_xq2_zpt2phi");
  if (!hGen) {
    std::cerr << "[ERROR] no TH2D found in gen file\n";
    return 4;
  }

  // load maps
  auto mapData = LoadFitMap(tData, "DATA (rec_data)");
  auto mapRec  = LoadFitMap(tRec,  "REC  (rec_true)");

  // create output
  TFile fout(outp.c_str(), "RECREATE");
  if (fout.IsZombie()) { std::cerr << "[ERROR] cannot create " << outp << "\n"; return 5; }

  // output hists
  auto* hScaled = dynamic_cast<TH2D*>(hGen->Clone("h2_gen_scaled"));
  hScaled->Reset();
  hScaled->SetTitle("gen * nPions(data)/nPions(rec);x-bin (ix);y-bin (iy)");

  auto* hScale = dynamic_cast<TH2D*>(hGen->Clone("h2_scale_data_over_rec"));
  hScale->Reset();
  hScale->SetTitle("scale = nPions(data)/nPions(rec);x-bin (ix);y-bin (iy)");

  // output tree
  TTree tout("scale_tree", "per-(ix,iy): gen * nPions_data / nPions_rec");
  int ix=0, iy=0;
  double gen=0;
  double nData=0, eData=0;
  double nRec=0,  eRec=0;
  double scale=0;
  double value=0;
  double err_value=0;
  int ok_data=0, ok_rec=0, ok_all=0;

  tout.Branch("ix", &ix, "ix/I");
  tout.Branch("iy", &iy, "iy/I");
  tout.Branch("gen", &gen, "gen/D");
  tout.Branch("nData", &nData, "nData/D");
  tout.Branch("eData", &eData, "eData/D");
  tout.Branch("nRec", &nRec, "nRec/D");
  tout.Branch("eRec", &eRec, "eRec/D");
  tout.Branch("scale", &scale, "scale/D");
  tout.Branch("value", &value, "value/D");
  tout.Branch("err_value", &err_value, "err_value/D");
  tout.Branch("ok_data", &ok_data, "ok_data/I");
  tout.Branch("ok_rec", &ok_rec, "ok_rec/I");
  tout.Branch("ok_all", &ok_all, "ok_all/I");

  // loop bins (use same bin indices as TH2)
  const int nx = hGen->GetNbinsX();
  const int ny = hGen->GetNbinsY();

  // if gen has Sumw2 it can provide errors; if not, errors are 0
  if (hGen->GetSumw2N() == 0) hGen->Sumw2();

  for (int bix = 1; bix <= nx; ++bix) {
    for (int biy = 1; biy <= ny; ++biy) {
      ix = bix;
      iy = biy;

      gen = hGen->GetBinContent(bix, biy);
      const double genErr = hGen->GetBinError(bix, biy);

      auto itD = mapData.find({bix,biy});
      auto itR = mapRec.find({bix,biy});

      ok_data = (itD != mapData.end()) ? 1 : 0;
      ok_rec  = (itR != mapRec.end())  ? 1 : 0;

      nData = ok_data ? itD->second.n   : -1.0;
      eData = ok_data ? itD->second.err : -1.0;
      nRec  = ok_rec  ? itR->second.n   : -1.0;
      eRec  = ok_rec  ? itR->second.err : -1.0;

      ok_all = (ok_data && ok_rec && (nRec > 0.0) && (nData > 0.0)) ? 1 : 0;

      if (!ok_all || gen == 0.0) {
        scale = 0.0;
        value = 0.0;
        err_value = 0.0;
        hScaled->SetBinContent(bix, biy, 0.0);
        hScale->SetBinContent(bix, biy, 0.0);
        tout.Fill();
        continue;
      }

      scale = nData / nRec;
      value = gen * scale;

      // optional propagated error (includes gen bin error if present)
      const double relD = (eData > 0.0 && nData > 0.0) ? (eData / nData) : 0.0;
      const double relR = (eRec  > 0.0 && nRec  > 0.0) ? (eRec  / nRec ) : 0.0;
      const double relG = (genErr > 0.0 && gen > 0.0) ? (genErr / gen) : 0.0;
      err_value = std::abs(value) * std::sqrt(relG*relG + relD*relD + relR*relR);

      hScaled->SetBinContent(bix, biy, value);
      hScaled->SetBinError(bix, biy, err_value);

      hScale->SetBinContent(bix, biy, scale);

      tout.Fill();
    }
  }

  fout.cd();
  hGen->Write("h2_gen_input"); // keep a copy for reference
  hScale->Write();
  hScaled->Write();
  tout.Write();
  fout.Write();
  fout.Close();

  std::cout << "[OK] wrote: " << outp << "\n";
  return 0;
}

// ROOT driver
void run_build_gen_scaled(const char* out_path = nullptr) {
  std::cout << "[RUN] data fit: " << kDataFitFile << "\n"
            << "[RUN] rec  fit: " << kRecFitFile  << "\n"
            << "[RUN] gen hist: " << kGenHistsFile << "\n";
  BuildGenScaled(out_path);
}

#ifndef __CLING__
int main(int argc, char** argv) {
  const char* outp = (argc >= 2) ? argv[1] : kDefaultOut;
  return BuildGenScaled(outp);
}
#endif
