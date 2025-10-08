// compile: root-config --cxx --cflags -std=c++17 rdf_all_in_dir.cxx $(root-config --libs) -o rdf_all_in_dir
#include <ROOT/RDataFrame.hxx>
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TString.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <string>
#include <ROOT/TThreadExecutor.hxx>   // for GetThreadPoolSize
#include <TGraphErrors.h>  // ADD THIS
#include <chrono>
#include <iomanip>

#include "binning_params.cxx" // assumes this defines: N_pi0mm_bins, peak_half_width, etc.
using namespace std;

const double peak_half_width = 0.06;
const double minPi0Mass = 0.133 - peak_half_width;
const double maxPi0Mass = 0.133 + peak_half_width;

double GetHistMaxInRange(TH1D *inhist,double minMass, double maxMass){
  const int sp = inhist->FindBin(minMass);
  double maxV = inhist->GetBinContent(sp);
  for (int i_bin = sp; i_bin <= inhist->FindBin(maxMass); i_bin++){
    double curV = inhist->GetBinContent(i_bin);
    if (curV > maxV) maxV = curV;
  }
  return maxV;
}

// ---- key for the 4D bin combination ----
struct Key {
  int xq2, xq2g, zpt2, zpt2g;
  bool operator==(const Key& o) const {
    return xq2==o.xq2 && xq2g==o.xq2g && zpt2==o.zpt2 && zpt2g==o.zpt2g;
  }
};

// hash for unordered_map
struct KeyHash {
  std::size_t operator()(const Key& k) const noexcept {
    std::size_t h = 1469598103934665603ull;
    auto mix = [&](uint64_t v){ h ^= v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2); };
    mix((uint64_t)k.xq2);  mix((uint64_t)k.xq2g);
    mix((uint64_t)k.zpt2); mix((uint64_t)k.zpt2g);
    return h;
  }
};

// make a unique hist name per slot+key
static std::string HistName(unsigned slot, const Key& k) {
  std::ostringstream os;
  os << "h_s" << slot << "_xq2_" << k.xq2 << "_xq2g_" << k.xq2g
     << "_zpt2_" << k.zpt2 << "_zpt2g_" << k.zpt2g;
  return os.str();
}



double GetLastNonTrashBinCenter(TH1D *inhist){
  const double start_mass = 0.13;
  double last_bin_content = inhist->GetBinContent(inhist->FindBin(start_mass));
  for (int iBin = inhist->FindBin(start_mass); iBin < inhist->GetNbinsX() - 2; iBin++){
    double cur_bin_con = inhist->GetBinContent(iBin);
    if (last_bin_content + cur_bin_con + inhist->GetBinContent(iBin + 1) < 6)
      return inhist->GetBinCenter(iBin + 1); 
    last_bin_content = cur_bin_con;
  }
  return inhist->GetBinCenter(inhist->GetNbinsX() - 1); 
}


static std::string fmt_hms(long long sec_total) {
  int h = static_cast<int>(sec_total / 3600);
  int m = static_cast<int>((sec_total % 3600) / 60);
  int s = static_cast<int>(sec_total % 60);
  std::ostringstream os;
  os << h << ':' << std::setw(2) << std::setfill('0') << m
     << ':' << std::setw(2) << std::setfill('0') << s;
  return os.str();
}

static void print_progress(std::size_t i, std::size_t total,
                           std::chrono::steady_clock::time_point t0) {
  using namespace std::chrono;
  if (total == 0) return;
  double frac = static_cast<double>(i) / static_cast<double>(total);
  int barw = 40;
  int filled = static_cast<int>(frac * barw);
  auto elapsed = duration_cast<seconds>(steady_clock::now() - t0).count();
  double rate = (elapsed > 0) ? (static_cast<double>(i) / elapsed) : 0.0;
  long long eta = (rate > 0.0) ? static_cast<long long>((total - i) / rate) : 0;

  std::ostringstream os;
  os << '\r' << '[';
  for (int k = 0; k < barw; ++k) os << (k < filled ? '=' : ' ');
  os << "] " << std::fixed << std::setprecision(1) << (frac * 100.0) << "%  "
     << i << '/' << total << "  "
     << "elapsed " << fmt_hms(elapsed) << "  "
     << "ETA " << fmt_hms(eta) << ' ';
  std::cout << os.str() << std::flush;
}



void fit_and_fill_response(const std::string path_to_folder_IN) {
  ROOT::EnableImplicitMT(); // comment out for single-thread

  const std::string dir  = path_to_folder_IN;
  const std::string tree = "h22";

  ROOT::RDataFrame df(tree, dir + "/*.root");

  std::cout << "Entries: " << *df.Count() << "\n";
  for (auto& c : df.GetColumnNames()) std::cout << "  " << c << "\n";

  auto d = df;

  // Per-thread (slot) maps: Key -> TH1D
  unsigned nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;

  using HMap = std::unordered_map<Key, std::unique_ptr<TH1D>, KeyHash>;
  std::vector<HMap> slot_maps(nSlots);

  // Fill per-slot histograms
  d.ForeachSlot(
    [&](unsigned slot, int xq2, int xq2g, int zpt2, int zpt2g, double MM) {
      Key key{xq2, xq2g, zpt2, zpt2g};
      auto& m = slot_maps[slot];
      auto it = m.find(key);
      if (it == m.end()) {
        auto hname  = HistName(slot, key);
        auto htitle = ";MM;entries";
        auto h = std::make_unique<TH1D>(hname.c_str(), htitle,
                                        N_pi0mm_bins, minPi0Mass, maxPi0Mass); // FIX names
        h->SetDirectory(nullptr);
        TH1::SetDefaultSumw2();
        h->Sumw2();
        h->SetCanExtend(TH1::kAllAxes);
        it = m.emplace(key, std::move(h)).first;
      }
      it->second->Fill(MM);
    },
    {"bin_xBQ2_Valerii","bin_xBQ2_Valeriim","zpt2phit_8x8x9","zpt2phit_8x8x9m","pi0_m"}
  );

  // Merge slot maps into one
  HMap merged;
  for (unsigned s = 0; s < nSlots; ++s) {
    for (auto& kv : slot_maps[s]) {
      const Key& key = kv.first;
      TH1D* hslot = kv.second.get();
      auto it = merged.find(key);
      if (it == merged.end()) {
        auto hclone = std::unique_ptr<TH1D>(
            static_cast<TH1D*>(hslot->Clone(HistName(9999, key).c_str())));
        hclone->SetDirectory(nullptr);
        merged.emplace(key, std::move(hclone));
      } else {
        it->second->Add(hslot);
      }
    }
  }

  // Output
  const char* outFile = "gaus_fits.root";
  TFile fout(outFile, "RECREATE");
  TTree tout("h22_fit", "Gaussian sigma per (xq2bin,xq2bin_gen,zpt2bin,zpt2bin_gen)");
  int    out_xq2, out_xq2g, out_zpt2, out_zpt2g;
  double out_sigma;
  tout.Branch("xq2bin",      &out_xq2);
  tout.Branch("xq2bin_gen",  &out_xq2g);
  tout.Branch("zpt2bin",     &out_zpt2);
  tout.Branch("zpt2bin_gen", &out_zpt2g);
  tout.Branch("sigma_gaus",  &out_sigma);

  TFile foutH("fit_hists.root", "RECREATE");
  
  gStyle->SetOptStat(0);
  gSystem->mkdir("fit_pngs", /*recursive=*/kTRUE);
  
  // --- progress setup ---
  const std::size_t total_groups = merged.size();
  std::size_t done = 0;
  auto t0 = std::chrono::steady_clock::now();
  auto last_print = t0;
  // ----------------------
  
  for (auto& kv : merged) {
    ++done;
  
    const Key& key = kv.first;
    TH1D* h = kv.second.get();

    double sigma = std::numeric_limits<double>::quiet_NaN();
    if (h->GetEntries() >= 30) {
      const double hardMaxFitEdge = 0.35;
      const double edgeOfevents   = GetLastNonTrashBinCenter(h);
      const double max_fit        = edgeOfevents > hardMaxFitEdge ? hardMaxFitEdge : edgeOfevents;

      // 1) Peak: Gaussian
      auto fcry = new TF1(("fcry_" + std::string(h->GetName())).c_str(),
                          "gaus", 0.07, edgeOfevents > 0.18 ? 0.18 : edgeOfevents);
      fcry->SetParameters(h->GetBinContent(h->FindBin(0.133)), 0.133, 0.011);
      fcry->SetParLimits(1, 0.133 - peak_half_width/2, 0.133 + peak_half_width/2);
      fcry->SetLineColor(kGreen+2);
      fcry->SetLineWidth(2);
      h->Fit(fcry, "QRS0"); // 0 = don't draw (we're saving, not displaying)

      // 2) Background: linear least squares on sidebands
      TGraphErrors g;
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double x  = h->GetBinCenter(i);
        const double y  = h->GetBinContent(i);
        const double ey = std::sqrt(std::max(0.0, y));
        if (x < minPi0Mass || x > maxPi0Mass) {
          const int n = g.GetN();
          g.SetPoint(n, x, y);
          g.SetPointError(n, 0.0, ey);
        }
      }
      auto bkg = new TF1(("bkg_" + std::string(h->GetName())).c_str(),
                         "pol3", 0.05, max_fit);            // HEAP object
      bkg->SetLineColor(kOrange+7);
      bkg->SetLineStyle(2);
      bkg->SetLineWidth(2);
      g.Fit(bkg, "QRS0");

      // 3) Combined: gaus(0) + pol3(3)
      auto fc_poly = new TF1(("fc_poly_" + std::string(h->GetName())).c_str(),
                             "gaus(0)+pol3(3)", 0.05, max_fit);
      double par[8] = {0};             // 3 (gaus) + 4 (pol3) = 7 params; keep a bit extra
      fcry->GetParameters(&par[0]);
      bkg->GetParameters(&par[3]);
      fc_poly->SetParameters(par);
      fc_poly->SetParLimits(1, 0.133 - peak_half_width, 0.133 + peak_half_width);
      fc_poly->SetNpx(500);
      fc_poly->SetLineColor(kRed);
      fc_poly->SetLineWidth(2);
      h->Fit(fc_poly, "QRS0");

      // overlay: attach functions to the histogram so they are saved with it
      h->GetListOfFunctions()->Add(fcry);
      h->GetListOfFunctions()->Add(bkg);
      h->GetListOfFunctions()->Add(fc_poly);

      // --- PNG export ---
      TCanvas c(Form("c_%s", h->GetName()), "", 900, 650);
      c.cd();
      
      // draw histogram with errors
      h->SetLineColor(kBlack);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(0.7);
      h->Draw("E");
      
      // draw the three fits on top
      fcry->Draw("same");
      bkg->Draw("same");
      fc_poly->Draw("same");
      
      // legend
      auto leg = new TLegend(0.55, 0.68, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(h,      "Data",            "lep");
      leg->AddEntry(fcry,   "Gaussian (peak)", "l");
      leg->AddEntry(bkg,    "Background (pol3)","l");
      leg->AddEntry(fc_poly,"Gaus+pol3 (sum)", "l");
      leg->Draw();
      
      // optional: logy if useful
      // c.SetLogy();
      
      c.SaveAs(Form("fit_pngs/%s.png", h->GetName()));
      // also save a PDF if you like:
      // c.SaveAs(Form("fit_pngs/%s.pdf", h->GetName()));
      // --- end PNG export ---


      

      sigma = fc_poly->GetParameter(2);   // Gaussian sigma
      using namespace std::chrono;
      auto now = steady_clock::now();
      if (done == total_groups || now - last_print >= 100ms) {
        print_progress(done, total_groups, t0);
        last_print = now;
      }
    }

    // Fill the TTree with per-group sigma
    out_xq2   = key.xq2;
    out_xq2g  = key.xq2g;
    out_zpt2  = key.zpt2;
    out_zpt2g = key.zpt2g;
    out_sigma = sigma;
    tout.Fill();

    // Write the histogram (with attached 3 fits) into the other ROOT file
    foutH.cd();
    h->Write();   // the attached TF1s go with the histogram
  }

  // finalize both outputs
  tout.Write();
  fout.Close();

  foutH.Write();
  foutH.Close();

  std::cout << "Wrote sigmas to " << outFile << " (TTree 'h22_fit')\n";
  std::cout << "Wrote histograms with 3 overlaid fits to fit_hists.root\n";
}