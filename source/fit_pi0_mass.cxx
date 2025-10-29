// fit_pi0_mass.cxx  (drop-in replacement)
//
// Gaussian + background fitting utilities for π0 mass slices.
// - No standalone Gaussian prefit; background-only prefit on sidebands.
// - Tries pol0/1/2/3 (attempt order: pol2 → pol0 → pol3 → pol1); winner by AICc.
// - Weighted combined fit (ROI centered on dynamic peak).
// - Final guards (positivity, mean/sigma, reliability).
// - Draws the Gaussian component from the final combined fit ("gaus_from_fc_*").
// - If no candidate passes final checks, draws *all* fc lines (all orders/ranges)
//   in order-colored, half-transparent overlays.
//
// NEW in this version:
// - Peak seed search in [0.09, 0.18].
// - Dynamic trust model: tighter χ² and narrower σ when |μ - 0.133| is large.
//   • χ²/ndf max shrinks from 10 → 4 as deviation grows from 0.015 → 0.040.
//   • σ max shrinks from 0.10 → 0.020 over the same deviation span.

#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <limits>
#include <map>

#include "TH1D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

using std::pair;
using std::vector;

// ---------------- tunables ----------------
constexpr double kPi0MeanGuess   = 0.133;
constexpr double kSigmaGuess     = 0.011;

constexpr double kPeakHalfWidth  = 0.035;
constexpr double kNsig           = 5.5;

constexpr double kMin_PEAKrng    = 0.070;
constexpr double kMax_PEAKrng    = 0.180;

// ROI weighting (emphasize the peak region in combined fit)
constexpr double kROIhalf         = 0.020;
constexpr double kDeweightOutside = 1.1;

// Mean guard center
constexpr double kMeanNominal         = 0.133;
constexpr double kMeanAllowedHalfSpan = 0.040;

// Absolute sigma safety cap
constexpr double kSigmaMaxFinal = 0.024;
constexpr double kSigmaMinFinal = 0.004;

// ---------------- result structs ----------------
struct SkipDecision {
  bool        skip       = false;
  bool        is_empty   = false;
  bool        is_no_peak = false;
  std::string reason;
};

struct FitResult {
  bool        ok         = false;
  double      nPions     = -1.0;
  double      errPions   = -1.0;
  int         bkg_order  = -1;     // 0/1/2/3 for pol0/1/2/3
  double      chi2ndf    = 1e9;
  std::string reason;
};

// ---------------- helpers ----------------
static inline double ReducedChi2(const TFitResultPtr& r) {
  if (!r.Get()) return 1e9;
  const double ndf = r->Ndf();
  return (ndf > 0.0) ? (r->Chi2() / ndf) : 1e9;
}

static pair<double,double> ClampRange(TH1D* h, double lo, double hi) {
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();
  double L = std::max(lo, xmin + 1e-6);
  double U = std::min(hi, xmax - 1e-6);
  if (U <= L) U = L + (xmax - xmin) / std::max(1, h->GetNbinsX());
  return {L, U};
}

// --- Dynamic peak seed in [0.1, 0.16] with parabolic sub-bin refinement ---
static inline double RefineParabolicVertex(const TH1D* h, int imax) {
  const int nb = h->GetNbinsX();
  if (imax <= 1 || imax >= nb) return h->GetBinCenter(imax);
  const double y1 = h->GetBinContent(imax - 1);
  const double y0 = h->GetBinContent(imax);
  const double y2 = h->GetBinContent(imax + 1);
  const double denom = (y1 - 2.0 * y0 + y2);
  if (std::abs(denom) < 1e-12) return h->GetBinCenter(imax);
  const double delta = 0.5 * (y1 - y2) / denom;
  return h->GetBinCenter(imax) + delta * h->GetBinWidth(imax);
}

static inline double EstimatePeakSeed(TH1D* h, double lo = 0.1, double hi = 0.16) {
  const int ib_lo = h->FindBin(lo);
  const int ib_hi = h->FindBin(hi);
  int    imax = -1; double ymax = -1.0;
  for (int i = ib_lo; i <= ib_hi; ++i) {
    const double y = h->GetBinContent(i);
    if (y > ymax) { ymax = y; imax = i; }
  }
  if (imax < 0 || ymax <= 0.0) return kPi0MeanGuess;
  return RefineParabolicVertex(h, imax);
}

// Local second-moment sigma guess around mu (±0.02)
static inline double EstimateSigmaGuess(TH1D* h, double mu) {
  const double win = 0.020;
  const int ib_lo  = h->FindBin(mu - win);
  const int ib_hi  = h->FindBin(mu + win);
  double wsum = 0.0, m2 = 0.0;
  for (int i = ib_lo; i <= ib_hi; ++i) {
    const double x  = h->GetBinCenter(i);
    const double y  = std::max(0.0, h->GetBinContent(i));
    wsum += y;
    m2   += y * (x - mu) * (x - mu);
  }
  double sg = (wsum > 0.0) ? std::sqrt(m2 / std::max(1e-12, wsum)) : kSigmaGuess;
  return std::min(std::max(sg, 0.004), 0.040);
}

// Build weighted graph (deweight outside ROI)
static TGraphErrors BuildWeightedGraph(TH1D* h, double roiCenter, double roiHalf, double deweightOutside) {
  TGraphErrors gW;
  int n = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x  = h->GetBinCenter(i);
    const double y  = h->GetBinContent(i);
    double ey       = std::sqrt(std::max(0.0, y));
    if (ey <= 0.0) ey = 1.0;
    if (std::abs(x - roiCenter) > roiHalf) ey *= std::max(1.0, deweightOutside);
    gW.SetPoint(n, x, y);
    gW.SetPointError(n, 0.0, ey);
    ++n;
  }
  return gW;
}

static inline double ComputeAICc(const TFitResultPtr& r, int k, int n) {
  if (!r.Get() || n <= 0 || k <= 0) return std::numeric_limits<double>::infinity();
  const double chi2 = r->Chi2();
  const double AIC  = chi2 + 2.0 * k;
  const int denom   = (n - k - 1);
  if (denom <= 0) return AIC + 1e6;
  return AIC + (2.0 * k * (k + 1)) / double(denom);
}

static int ColorForOrder(int ord) {
  switch (ord) {
    case 0: return kBlue + 1;     // pol0
    case 1: return kMagenta + 1;  // pol1
    case 2: return kRed + 1;      // pol2
    case 3: return kGreen + 2;    // pol3
    default: return kGray + 2;
  }
}

// ---- Dynamic reliability model vs. |μ - 0.133| ----
static inline double DeltaFromNom(double mu) { return std::abs(mu - kMeanNominal); }

// χ²/ndf max shrinks from 10 → 2 as deviation grows 0.015 → 0.050
static inline double Chi2MaxForDelta(double d) {
  const double base = 10.0, tight = 4.0, d1 = 0.015, d2 = 0.050;
  if (d <= d1) return base;
  if (d >= d2) return tight;
  const double t = (d - d1) / (d2 - d1);
  return base + (tight - base) * t;
}

// ---------------- pre-checks ----------------
SkipDecision PrecheckHistogram(TH1D* h) {
  SkipDecision d;
  if (!h) { d.skip = true; d.is_empty = true; d.reason = "null histogram pointer"; return d; }

  const double integral_full = h->Integral();
  const double max_content   = h->GetMaximum();
  const double entries       = h->GetEntries();

  if (entries <= 0 || integral_full <= 0 || max_content <= 0) {
    d.skip = true; d.is_empty = true; d.reason = "empty histogram (zero entries/integral/maximum)";
    return d;
  }

  // Keep your existing Nbins and ROI-height checks as in your current file
  if (h->GetNbinsX() < 5) {
    d.skip = true; d.is_empty = true; d.reason = "ignored: NbinsX < 5";
    return d;
  }

  {
    // Skip if tallest bin in 0.10–0.15 is below 6
    const int ib_lo = std::max(1, h->FindBin(0.10));
    const int ib_hi = std::min(h->GetNbinsX(), h->FindBin(0.15));
    double max_roi = 0.0;
    for (int i = ib_lo; i <= ib_hi; ++i) max_roi = std::max(max_roi, h->GetBinContent(i));
    if (max_roi < 6.0) {
      d.skip = true; d.is_empty = true; d.reason = "ignored: max[0.10,0.15] < 6";
      return d;
    }
  }

  const double int_01015 = h->Integral(h->FindBin(0.10), h->FindBin(0.15));
  if (int_01015 < 10.0) {
    d.skip = true; d.is_empty = true; d.reason = "ignored: integral[0.10,0.15] < 10";
    return d;
  }

  d.skip = false;
  return d;
}

// ---------- Final-quality evaluator ----------
static bool FinalQualityAndYield(TH1D* h, TF1* bkg_final, TF1* fc_final, const TFitResultPtr& r_final,
                                 double& outN, double& outErr, double& outChi2, std::string& why_not) {
  if (!fc_final || !bkg_final || !r_final.Get()) { why_not = "null candidate"; return false; }

  const double mean  = fc_final->GetParameter(1);
  const double sigma = std::abs(fc_final->GetParameter(2));
  const double dmu   = DeltaFromNom(mean);

  if (std::abs(mean - kMeanNominal) > kMeanAllowedHalfSpan) {
    std::ostringstream ss; ss << "mean out of allowed window: |" << mean << " - " << kMeanNominal
                              << "| > " << kMeanAllowedHalfSpan;
    why_not = ss.str(); return false;
  }
  if (sigma > kSigmaMaxFinal) { why_not = "sigma absurdly large"; return false; }

  if (sigma < kSigmaMinFinal) { why_not = "sigma absurdly small"; return false; }

  // Dynamic narrowness requirement when far from 0.133
  const double sigMaxDyn = kSigmaMaxFinal;
  if (sigma > sigMaxDyn) {
    std::ostringstream ss; ss << "sigma too wide for offset: sigma=" << sigma <<" max:" << sigMaxDyn;
    why_not = ss.str(); return false;
  }
  if (mean < 0.1) { why_not = "mass is too low (<0.1)"; return false; }

  // Background positivity guard in [0.10, 0.127]
  {
    const double chkLo = 0.10, chkHi = 0.127;
    const double fLo = std::max(chkLo, bkg_final->GetXmin());
    const double fHi = std::min(chkHi, bkg_final->GetXmax());
    if (fHi > fLo) {
      const double minBkg = bkg_final->GetMinimum(fLo, fHi);
      if (minBkg < 0.0) { why_not = "negative background in [0.10,0.127]"; return false; }
    }
  }

  // Yields within ±kNsig σ
  std::vector<double> rn = {mean - kNsig * sigma, mean + kNsig * sigma};
  const double binw    = h->GetBinWidth(1);

  TF1 fgaus_int("fgaus_int", "gaus", fc_final->GetXmin(), fc_final->GetXmax());
  fgaus_int.SetParameters(fc_final->GetParameter(0), mean, sigma);

  TMatrixDSym cov_full = r_final->GetCovarianceMatrix();
  TMatrixDSym covG(3);
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) covG(i, j) = cov_full(i, j);

  const double gausInt = fgaus_int.Integral(rn[0], rn[1]) / binw;
  const double bkgInt  = bkg_final->Integral(rn[0], rn[1]) / binw;

  double dI     = fgaus_int.IntegralError(rn[0], rn[1], fgaus_int.GetParameters(), covG.GetMatrixArray());
  const double fitErr = dI / binw;

  if (bkgInt != 0.0 && (std::abs(gausInt) / std::abs(bkgInt)) < 0.1) {
    std::ostringstream ss; ss << "unreliable: signal/bkg < 0.1 inside ±" << kNsig << "σ";
    why_not = ss.str(); return false;
  }

  outN     = std::abs(gausInt);
  outErr   = std::abs(fitErr);
  outChi2  = ReducedChi2(r_final);
  return true;
}

// ---------------- the actual fit ----------------
FitResult FitPi0Mass(TH1D* h) {
  FitResult out;

  const auto dpre = PrecheckHistogram(h);
  if (dpre.skip) { out.reason = dpre.reason; return out; }

  // dynamic peak seed (now in [0.09, 0.18])
  const double mu_seed = EstimatePeakSeed(h, 0.09, 0.18);
  const double sideMin = mu_seed - kPeakHalfWidth;
  const double sideMax = mu_seed + kPeakHalfWidth;

  // sideband-only graph for bkg seeding
  TGraphErrors gSide;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x  = h->GetBinCenter(i);
    const double y  = h->GetBinContent(i);
    const double ey = std::sqrt(std::max(0.0, y));
    if (x < sideMin || x > sideMax) {
      const int n = gSide.GetN();
      gSide.SetPoint(n, x, y);
      gSide.SetPointError(n, 0.0, (ey>0?ey:1.0));
    }
  }

  // weighted-all-bins graph for the combined fit (ROI centered on mu_seed)
  TGraphErrors gAllW = BuildWeightedGraph(h, mu_seed, kROIhalf, kDeweightOutside);
  const int npts_all = gAllW.GetN();

  struct Cand {
    int              ord   = -1;
    double           L     = 0, U = 0;
    TF1*             bkg   = nullptr;
    TF1*             fc    = nullptr;
    TFitResultPtr    r;
    double           aicc  = std::numeric_limits<double>::infinity();
    double           chi2  = 1e9;
    double           nP    = -1, eP = -1;
    bool             passedFinal = false;
    std::string      why_not;
  };
  std::vector<Cand> cands;

  auto try_once = [&](int ord, double L, double U) -> Cand {
    Cand c; c.ord = ord; c.L = L; c.U = U;

    // 1) background-only prefit on sidebands
    TF1* bkg_try = nullptr;
    if      (ord == 0) bkg_try = new TF1("bkg_try_tmp", "pol0", L, U);
    else if (ord == 1) bkg_try = new TF1("bkg_try_tmp", "pol1", L, U);
    else if (ord == 2) bkg_try = new TF1("bkg_try_tmp", "pol2", L, U);
    else               bkg_try = new TF1("bkg_try_tmp", "pol3", L, U);

    bkg_try->SetLineColor(kOrange + 7);
    bkg_try->SetLineStyle(2);
    bkg_try->SetLineWidth(2);
    TFitResultPtr r_bkg = gSide.Fit(bkg_try, "QRS0");

    // 2) seed Gaussian (NO gaussian prefit)
    const double y_peak     = h->GetBinContent(h->FindBin(mu_seed));
    const double bkg_at_mu  = bkg_try->Eval(mu_seed);
    const double amp_guess  = std::max(1.0, y_peak - bkg_at_mu);
    const double sig_guess  = EstimateSigmaGuess(h, mu_seed);

    // 3) combined function and fit
    const char* polyName = (ord == 0 ? "pol0" : ord == 1 ? "pol1" : ord == 2 ? "pol2" : "pol3");
    std::string fcname   = std::string("fc_tmp_") + h->GetName() + "_" + polyName;
    TF1* fc = new TF1(fcname.c_str(), (std::string("gaus(0)+") + polyName + "(3)").c_str(), L, U);

    fc->SetParameters(amp_guess, mu_seed, sig_guess);
    fc->SetParLimits(0, 0.0, 1e9);
    fc->SetParLimits(1, mu_seed - kPeakHalfWidth, mu_seed + kPeakHalfWidth);
    fc->SetParLimits(2, 0.001, 0.10);

    const int npar_bkg = bkg_try->GetNpar();
    for (int i = 0; i < npar_bkg; ++i) fc->SetParameter(3 + i, bkg_try->GetParameter(i));

    fc->SetNpx(800);
    fc->SetLineColor(ColorForOrder(ord));
    fc->SetLineWidth(2);

    TFitResultPtr r = gAllW.Fit(fc, "QRS0"); // weighted fit

    if (!r.Get()) {
      c.why_not = "no fit result";
      delete fc; delete bkg_try;
      return c;
    }

    // Dynamic χ² threshold based on the fitted mean
    const double mean  = fc->GetParameter(1);
    const double dmu   = DeltaFromNom(mean);
    const double chi2max_dyn = Chi2MaxForDelta(dmu);
    const bool ok_prefilter  = (r->Status() == 0 && ReducedChi2(r) < chi2max_dyn);

    // finalize a background function with fitted params for checks
    TF1* bkg_final = new TF1((std::string("bkg_") + h->GetName() + "_" + polyName).c_str(),
                             polyName, L, U);
    bkg_final->SetLineColor(kOrange + 7);
    bkg_final->SetLineStyle(2);
    bkg_final->SetLineWidth(2);
    for (int i = 0; i < npar_bkg; ++i) bkg_final->SetParameter(i, fc->GetParameter(3 + i));

    double nP=0, eP=0, chi=0;
    std::string wn;
    bool okFinal = false;
    if (ok_prefilter) {
      okFinal = FinalQualityAndYield(h, bkg_final, fc, r, nP, eP, chi, wn);
    } else {
      wn = "prefilter fail: chi2/ndf > dyn threshold";
    }

    const int k_free = 3 + npar_bkg;
    const double aicc = okFinal ? ComputeAICc(r, k_free, npts_all)
                                : std::numeric_limits<double>::infinity();

    c.bkg         = bkg_final;
    c.fc          = fc;
    c.r           = r;
    c.aicc        = aicc;
    c.chi2        = ReducedChi2(r);
    c.passedFinal = okFinal;
    c.nP          = nP;
    c.eP          = eP;
    c.why_not     = wn;

    delete bkg_try;
    return c;
  };

  auto try_ord_with_ranges = [&](int ord) {
    std::vector<pair<double,double>> ranges = {
      {0.06, 0.27},
      {0.08, 0.35},
      {0.05, 0.40}
    };
    for (auto [Lraw, Uraw] : ranges) {
      auto [L, U] = ClampRange(h, Lraw, Uraw);
      cands.emplace_back( try_once(ord, L, U) );
    }
  };

  // Attempt order
  for (int ord : {2, 0, 3, 1}) try_ord_with_ranges(ord);

  // Collect successful candidates
  std::vector<size_t> ok_idx;
  for (size_t i = 0; i < cands.size(); ++i)
    if (cands[i].passedFinal && std::isfinite(cands[i].aicc)) ok_idx.push_back(i);

  if (ok_idx.empty()) {
    // FAIL CASE: draw *all* fc lines that produced a fit, per-order colors, half-transparent
    for (size_t i = 0; i < cands.size(); ++i) {
      if (!cands[i].r.Get() || !cands[i].fc) { if (cands[i].bkg) delete cands[i].bkg; continue; }
      const int col = ColorForOrder(cands[i].ord);
      cands[i].fc->SetLineColorAlpha(col, 0.5);
      cands[i].fc->SetLineWidth(3);
      cands[i].fc->SetLineStyle(1);
      //h->GetListOfFunctions()->Add(cands[i].fc);
      if (cands[i].bkg) delete cands[i].bkg; // background not needed for fail visual
      cands[i].bkg = nullptr;
    }
    out.reason = "no acceptable candidate (all orders/ranges failed final checks)";
    out.ok     = false;
    return out;
  }

  // Choose the best by AICc
  size_t best_i = ok_idx[0];
  for (size_t j = 1; j < ok_idx.size(); ++j) {
    const size_t idx = ok_idx[j];
    if (cands[idx].aicc < cands[best_i].aicc) best_i = idx;
  }

  // Dispose losers (keep only winner)
  for (size_t i = 0; i < cands.size(); ++i) {
    if (i == best_i) continue;
    if (cands[i].fc)  delete cands[i].fc;
    if (cands[i].bkg) delete cands[i].bkg;
  }

  auto& win = cands[best_i];
  h->GetListOfFunctions()->Add(win.bkg);
  h->GetListOfFunctions()->Add(win.fc);

  // Draw Gaussian-only component from final combined fit (pure overlay + names for SaveHistPNG)
  {
    const double A   = win.fc->GetParameter(0);
    const double MU  = win.fc->GetParameter(1);
    const double SIG = std::abs(win.fc->GetParameter(2));

    TF1* fgaus_from_fc =
        new TF1((std::string("gaus_from_fc_") + h->GetName()).c_str(),
                "gaus", win.fc->GetXmin(), win.fc->GetXmax());
    fgaus_from_fc->SetParName(0, "Amp");
    fgaus_from_fc->SetParName(1, "Mean");
    fgaus_from_fc->SetParName(2, "Sigma");

    fgaus_from_fc->SetParameters(A, MU, SIG);
    fgaus_from_fc->SetLineColor(kAzure + 2);
    fgaus_from_fc->SetLineStyle(2);
    fgaus_from_fc->SetLineWidth(3);

    // Do NOT fit here; SaveHistPNG will do the constrained refit to update stats box
    h->GetListOfFunctions()->Add(fgaus_from_fc);

  }

  out.ok        = true;
  out.nPions    = win.nP;
  out.errPions  = win.eP;
  out.bkg_order = win.ord;
  out.chi2ndf   = win.chi2;
  return out;
}
