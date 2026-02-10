// fit_pi0_mass_valerii_adapter.cxx
//
// Adapter layer so split_and_fit_unified.cxx can keep using:
//   SkipDecision PrecheckHistogram(TH1D*)
//   FitResult    FitPi0Mass(TH1D*)
//
// Internally calls your Valerii fit:
//   fTuple findBestBack_Valerii(TH1D*)
//   int    badFit(fTuple, TH1D*)

#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
using namespace std; // or std::tuple/std::vector/std::get

#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"

// ---------------- result structs (match your unified driver) ----------------
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
  int         bkg_order  = -1;     // 0/1/2/3... (pol0..polN)
  double      chi2ndf    = 1e9;
  std::string reason;
};

// ---------------- helpers ----------------
static inline double ReducedChi2(const TFitResultPtr& r) {
  if (!r.Get()) return 1e9;
  const double ndf = r->Ndf();
  return (ndf > 0.0) ? (r->Chi2() / ndf) : 1e9;
}

// Keep consistent with your old driver’s behavior:
SkipDecision PrecheckHistogram(TH1D* h) {
  SkipDecision d;
  if (!h) { d.skip = true; d.is_empty = true; d.reason = "null histogram"; return d; }

  const double integral_full = h->Integral();
  const double max_content   = h->GetMaximum();
  const double entries       = h->GetEntries();

  if (entries <= 0 || integral_full <= 0 || max_content <= 0) {
    d.skip = true; d.is_empty = true;
    d.reason = "empty histogram (zero entries/integral/maximum)";
    return d;
  }

  if (h->GetNbinsX() < 5) {
    d.skip = true; d.is_empty = true;
    d.reason = "ignored: NbinsX < 5";
    return d;
  }

  // “No peak / too low stats” gate similar to what you already used
  {
    const int ib_lo = std::max(1, h->FindBin(0.10));
    const int ib_hi = std::min(h->GetNbinsX(), h->FindBin(0.15));
    double max_roi = 0.0;
    for (int i = ib_lo; i <= ib_hi; ++i) max_roi = std::max(max_roi, h->GetBinContent(i));
    if (max_roi < 6.0) {
      d.skip = true; d.is_no_peak = true;
      d.reason = "ignored: max[0.10,0.15] < 6";
      return d;
    }
  }

  const double int_01015 = h->Integral(h->FindBin(0.10), h->FindBin(0.15));
  if (int_01015 < 10.0) {
    d.skip = true; d.is_empty = false;
    d.reason = "ignored: integral[0.10,0.15] < 10";
    return d;
  }

  return d;
}

// ---------------- the actual fit wrapper ----------------
FitResult FitPi0Mass(TH1D* h) {
  FitResult out;

  const auto dpre = PrecheckHistogram(h);
  if (dpre.skip) { out.reason = dpre.reason; return out; }

  // Run Valerii fit (your code)
  fTuple fvec = findBestBack_Valerii(h);
  TF1* fSig   = std::get<0>(fvec);
  TF1* fBkg   = std::get<1>(fvec);
  TF1* fComb  = std::get<2>(fvec);
  TFitResultPtr r = std::get<3>(fvec);

  if (!fSig || !fBkg || !fComb || !r.Get()) {
    out.ok = false;
    out.reason = "Valerii fit returned null function or null fit result";
    return out;
  }

  // Decide OK/FAIL using your badFit() logic
  const int bf = badFit(fvec, h);
  const std::string hname = h->GetName();

  // IMPORTANT: give per-hist unique names so SaveHistPNG finds them
  // Also: keep only the final 3 functions on the histogram
  if (auto* lf = h->GetListOfFunctions()) lf->Clear("nodelete");

  TF1* gaus = (TF1*)fSig->Clone((std::string("gaus_from_fc_") + hname).c_str());
  TF1* bkg  = (TF1*)fBkg->Clone((std::string("bkg_")         + hname).c_str());
  TF1* fc   = (TF1*)fComb->Clone((std::string("fc_")          + hname).c_str());

  // Cosmetic defaults (optional)
  gaus->SetLineStyle(2);
  gaus->SetLineWidth(3);

  bkg->SetLineStyle(2);
  bkg->SetLineWidth(2);

  fc->SetLineWidth(2);

  h->GetListOfFunctions()->Add(bkg);
  h->GetListOfFunctions()->Add(fc);
  h->GetListOfFunctions()->Add(gaus);

  // polynomial order = (Npar(bkg) - 1) for polN
  out.bkg_order = std::max(-1, bkg->GetNpar() - 1);
  out.chi2ndf    = ReducedChi2(r);

  if (bf != 0) {
    std::ostringstream ss;
    ss << "badFit=" << bf;
    out.ok = false;
    out.reason = ss.str();
    return out;
  }

  // Yield + error from gaussian integral in ±Nσ (match your old driver behavior)
  constexpr double kNsig = 5.5;

  const double A    = fc->GetParameter(0);
  const double mean = fc->GetParameter(1);
  const double sig  = std::abs(fc->GetParameter(2));
  if (!(std::isfinite(mean) && std::isfinite(sig)) || sig <= 0.0) {
    out.ok = false;
    out.reason = "non-finite mean/sigma from fit";
    return out;
  }

  const double lo = mean - kNsig * sig;
  const double hi = mean + kNsig * sig;
  const double binw = h->GetBinWidth(1);

  TF1 fgaus_int("fgaus_int", "gaus", fc->GetXmin(), fc->GetXmax());
  fgaus_int.SetParameters(A, mean, sig);

  // covariance: take 3x3 gaussian block from full covariance
  TMatrixDSym cov_full = r->GetCovarianceMatrix();
  TMatrixDSym covG(3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      covG(i, j) = cov_full(i, j);

  const double I = fgaus_int.Integral(lo, hi) / binw;
  const double dI = fgaus_int.IntegralError(
      lo, hi,
      fgaus_int.GetParameters(),
      covG.GetMatrixArray()
  ) / binw;

  out.ok       = true;
  out.nPions   = std::abs(I);
  out.errPions = std::abs(dI);
  out.reason   = "ok";
  return out;
}


//Function tat sets histogram plot colors and points
void setPlotPar(TH1D &hist, Color_t col = kBlue){
     hist.SetLineColor(col);
     hist.SetLineWidth(2);
     hist.SetMarkerStyle(20);
     hist.SetMarkerColor(col);
}

//Type define for fitting, its used to return fit functions and fit parameters
typedef tuple<TF1*, TF1*, TF1*, TFitResultPtr> fTuple;

//Function fits a mass distribution and returns the integral scaled by 1/binwidth
//The fit uses a gaussian for the signal and a pol4 for the background
double fitIntegral(TH1D *inhist){

    double par[20];

    auto fcry    = new TF1("fcry", "gaus", 0.08, 0.35);//0.08-0.35
    auto fpoly   = new TF1("fpoly", "pol4", 0.08, 0.35);
    auto fc_poly = new TF1("fc_poly","gaus(0) +pol4(3)", 0.08, 0.35);

    fcry->SetParameters(1e5, 0.133, 0.011);
    fcry->SetParLimits(1, 0.1, 0.15);
    fc_poly->SetParLimits(1, 0.13, 0.14);
    fc_poly->SetParLimits(2, 0.01, 0.014);

    inhist->Fit(fcry, "QR");
    inhist->Fit(fpoly, "QR");
    fcry->GetParameters(&par[0]);

    fpoly->GetParameters(&par[3]);
    fc_poly->SetParameters(par);
    inhist->Fit(fc_poly, "QR");

    double pf[8];
    for (int j = 0; j < 8; j++){pf[j] = fc_poly->GetParameter(j);}
    fcry->SetParameters(&pf[0]);
    fpoly->SetParameters(&pf[3]);

    double mean  = fc_poly->GetParameter(1);
    double sigma = fc_poly->GetParameter(2);
    vector<double> range = {mean - n_sigma_pi0mass*sigma, mean + n_sigma_pi0mass*sigma};

    return fcry->Integral(range[0], range[1]) / inhist->GetBinWidth(1);
}

//Function fits a mass distribution and returns the rsulting fit functions and fit pointer
//The fit uses a gaussian for the signal and a pol4 for the background as default
//The pol degree can be changed by altering the pint variable;
fTuple fitFunc(TH1D *inhist, int pint = 4){

    double par[20];
    TString bstr = Form("pol%d", pint);
    TString cstr = Form("gaus(0) + %s(3)", bstr.Data());
    auto fcry    = new TF1("fcry",    "gaus",      0.08, 0.35);
    auto fpoly   = new TF1("fpoly",   bstr.Data(), 0.08, 0.35);
    auto fc_poly = new TF1("fc_poly", cstr.Data(), 0.08, 0.35);

    fcry->SetParameters(inhist->GetMaximum(), 0.133, 0.011);
    // peak position limits
    // gaus only
    fcry->SetParLimits(1, 0.1, 0.15);
    // gaus + pol: peak position
    fc_poly->SetParLimits(1, 0.12, 0.14);
    //gaus + pol: width
    fc_poly->SetParLimits(2, 0.08, 0.014);

    // fit gaus and pol separately
    auto r1 = inhist->Fit(fcry, "QRS");
    auto r2 = inhist->Fit(fpoly, "QRS");
    fcry->GetParameters(&par[0]);

    fpoly->GetParameters(&par[3]);
    // fit parametrs are set from gaus and pol 4 fits
    fc_poly->SetParameters(par);
    TFitResultPtr r3 = inhist->Fit(fc_poly, "QRS");

    // contains all the final fit parameters
    double pf[8];
    for (int j = 0; j < 8; j++){pf[j] = fc_poly->GetParameter(j);}
    fcry->SetParameters(&pf[0]);
    fpoly->SetParameters(&pf[3]);

    fcry->SetLineColor(kGreen);
    fpoly->SetLineColor(kOrange);
    fc_poly->SetLineColor(kRed);
    
    return (fTuple){fcry, fpoly, fc_poly, r3};
}

//Function fits a phi trento distribution with two functions:
//A(1+Bcos(phi) + Ccos(2phi))
//pol0
//The two functions are then returned in a vector
vector<TF1*> fitTrento(TH1D *inhist){
    double par[20];
    auto fn = new TF1("f", "[0]*(1 + [1]*cos(x*TMath::Pi()/180.0) + [2]*cos(2*x*TMath::Pi()/180.0))", 0, 360);//0.08-0.25
    fn->SetParLimits(0, 0, inhist->Integral());
    int nbins = inhist->GetNbinsX();

    double min = 0, max = 360;
    inhist->Fit(fn, "QR", "", min, max);//+rob=0.8
    fn->SetLineColor(kRed);

    auto fn2 = new TF1("l", "[0]", 0, 360);
    double maxl = inhist->GetBinContent(inhist->GetMaximumBin());
    double minl = inhist->GetBinContent(inhist->GetMinimumBin());
    double low = 0.2*maxl ;
    fn2->SetParameter(0, inhist->GetMean());
    //fn2->SetParLimits(0, low, maxl);
    inhist->Fit(fn2, "QR");
    fn2->SetLineColor(kGreen);

    return (vector<TF1*>){fn, fn2};
}

//Function fits a phi trento distribution and returns the A term and its associated error from
//A(1+Bcos(phi) + Ccos(2phi))
vector<double> fitTrentoA(TH1D *inhist){

    double par[20];
    auto fn = new TF1("f", "[0]*(1 + [1]*cos(x*TMath::Pi()/180.0) + [2]*cos(2*x*TMath::Pi()/180.0))", 0, 360);

    inhist->Fit(fn, "QR");
   
    return (vector<double>){fn->GetParameter(0), fn->GetParError(0)};
}

//Function fits a phi trento distribution with a pol0 and returns the A term and its associated error
vector<double> fitLinear(TH1D *inhist){

    double par[20];
    auto fn = new TF1("l", "pol0", 0, 360);

    inhist->Fit(fn, "QR");
   
    return (vector<double>){fn->GetParameter(0), fn->GetParError(0)};
}

//Functions returns a number is a mass fit is off
//if everything okay it returs 0 else it returns numbers
//related to what issue it is having
int badFit(fTuple fv, TH1D* h){
  int res = 0;
  double min = 0.09, max = 0.16;//changed from 0.1-0.168
  double h_int = h->Integral(h->FindBin(min), h->FindBin(max));
  double binw  = h->GetBinWidth(1);
  double f_int = get<2>(fv)->Integral(min, max) / binw;
  double s_int = get<0>(fv)->Integral(min, max) / binw;
  double b_int = get<1>(fv)->Integral(min, max) / binw;   
  if (h_int <= 0)   {res+=1;}    //number of events within 0.08-0.2 is <= 0
  if (s_int <= 0)   {res+=10;}   //No signal
  if (s_int > f_int){res+=100;}  //Signal > entire integral
  if (b_int < 0)    {res+=1000;} //Background is negative in signal range
  // Valerii
  if (h->Integral(h->FindFixBin(0.09), h->FindFixBin(0.15)) < 40) {res+=10000;} // statistic is too low
  if (get<0>(fv)->GetParameter(2) > 0.1) {res+=100000;} // peak is too wides
  return res;
}

//Function reduces the degree of the polynomial if there is a fitting 
//issue related to badFit
fTuple findBestBack(TH1D* inhist){
  auto fvec = fitFunc(inhist);
  TString istr = Form("%d", badFit(fvec, inhist));
  // no events, no way to improve, proceed
  if (istr.EndsWith("1")){return fvec;}
  if (istr == "0"){return fvec;}
  //if (istr == "10"){return fvec;}//Added To reflect above, s_int change
  // trying various pols
  for (int i = 3; i > -1; i--){
    fvec = fitFunc(inhist, i); 
    istr = Form("%d", badFit(fvec, inhist));
    if (istr.Atoi() == 0){break;}
  }
  if (istr.Atoi() == 0){return fvec;}
  return fvec;
}

int GetBkgrN(TH1D *inhist){
  const double cutOff_for_pol2 = 0.6;
  const double cutOff_for_pol1 = 0.8;
  const double minMass = 0.09, maxMass = 0.16;
  const double minFullRange = 0.05, maxFullRange = 0.35;
  
  double peakRegionInt = inhist->Integral(inhist->FindBin(minMass), inhist->FindBin(maxMass));
  double wholeRangeInt = inhist->Integral(inhist->FindBin(minFullRange), inhist->FindBin(maxFullRange));
  const double ratio = (peakRegionInt / wholeRangeInt);
  if (ratio > cutOff_for_pol1) return 0;
  if (ratio > cutOff_for_pol2) return 1;
  return 2;
}

double GetLastNonTrashBinCenter(TH1D *inhist){
  const double start_mass = 0.13;
  double last_bin_content = inhist->GetBinContent(inhist->FindBin(start_mass));
  for (int iBin = inhist->FindBin(start_mass); iBin < inhist->GetNbinsX() - 2; iBin++){
    double cur_bin_con = inhist->GetBinContent(iBin);
    if (last_bin_content + cur_bin_con + inhist->GetBinContent(iBin + 1) < 6) return inhist->GetBinCenter(iBin + 1); 
    last_bin_content = cur_bin_con;
    }
  return inhist->GetBinCenter(inhist->GetNbinsX() - 1); 
}

double GetHistMaxInRange(TH1D *inhist,double minMass, double maxMass){
  const int sp = inhist->FindBin(minMass);
  double maxV =  inhist->GetBinContent(sp);
  for (int i_bin = sp; i_bin <= inhist->FindBin(maxMass); i_bin++){
    double curV = inhist->GetBinContent(i_bin);
    if ( curV > maxV) maxV = curV;
  }
  return maxV;
}

// Valerii's fit:
fTuple fitFunc_Valerii(TH1D *inhist, const int pint = 2){

    // is used to delete the peak for bckrd fit
    // and for gaus peak param limits
    const double peak_half_width = 0.06;
  
    double par[20];
    TString bstr = Form("pol%d", pint);
    TString cstr = Form("gaus(0) + %s(3)", bstr.Data());

    const double hardMaxFitEdge = 0.35;
    const double edgeOfevents = GetLastNonTrashBinCenter(inhist);
    const double max_fit = edgeOfevents > hardMaxFitEdge ? hardMaxFitEdge : edgeOfevents;
  
    auto fcry    = new TF1("fcry",    "gaus",      0.07, edgeOfevents > 0.18 ? 0.18 : edgeOfevents);
    auto fpoly   = new TF1("fpoly",   bstr.Data(), 0.05, max_fit);
    auto fc_poly = new TF1("fc_poly", cstr.Data(), 0.05, max_fit);

    // prepare gaus fit
    fcry->SetParameters(inhist->GetBinContent(GetHistMaxInRange(inhist, 0.133 - peak_half_width, 0.133 + peak_half_width)), 0.133, 0.011);
    //fix peak position at correct value
    //fcry->FixParameter(1, 0.133);
    fcry->SetParLimits(1, 0.133 - peak_half_width/2, 0.133 + peak_half_width/2);
  
    // fit gaus and pol separately
    auto r1 = inhist->Fit(fcry, "QRS");
  
    // Fitting bckgr on a histogram were the peak is removed completely
    auto inHist_noPeak = (TH1D*)inhist->Clone();
    for (int iBin = inHist_noPeak->FindBin(0.133 - peak_half_width); 
          iBin <= inHist_noPeak->FindBin(0.133 + peak_half_width); iBin++){
      inHist_noPeak->SetBinContent(iBin, 0);
      inHist_noPeak->SetBinError(iBin, 0);
    } 
    auto r2 = inHist_noPeak->Fit(fpoly, "QRS");
  
    fcry->GetParameters(&par[0]);
    fpoly->GetParameters(&par[3]);
     
    // fit parametrs are set from gaus and pol N fits
    fc_poly->SetParameters(par);
    //fix peak position if it is out of the range
    //fc_poly->FixParameter(1, 0.133);
    fc_poly->SetParLimits(1, 0.133 - peak_half_width, 0.133 + peak_half_width);

    TFitResultPtr r3 = inhist->Fit(fc_poly, "QWWRS");

    // contains all the final fit parameters
    // gaus + polN
    const size_t nParams_full_fit_fun = 3 + (pint + 1);
    std::vector<double> pf(nParams_full_fit_fun);
    
    for (size_t j = 0; j < nParams_full_fit_fun; ++j)
      pf[j] = fc_poly->GetParameter(j);
    
    fcry->SetParameters(pf.data());       // params 0..2
    fpoly->SetParameters(pf.data() + 3);  // params 3..end

    fcry->SetLineColor(kGreen);
    fpoly->SetLineColor(kOrange);
    fc_poly->SetLineColor(kRed);
    
    return (fTuple){fcry, fpoly, fc_poly, r3};
}



fTuple findBestBack_Valerii(TH1D* inhist){
  int n_pol_bckgr = GetBkgrN(inhist);
  auto fvec = fitFunc_Valerii(inhist, n_pol_bckgr);
  TString istr = Form("%d", badFit(fvec, inhist));
  // no events, no way to improve, proceed
  if (istr.EndsWith("1")){return fvec;}
  // All good
  if (istr == "0"){return fvec;}
  
  // there is a problem
  // but still return the fit
  return fvec;
}


//Function takes in a 2D histogram with phi trento in the y axis.
//Fits the phi trento distribution in each x bin and returns a 
//1D histogram with the original x axis and the bin content and error
//being the extracted A and error on A from the phi trento fit.
TH1D* extractTrentoA(TH2D* inHist){
  int nbinsx = inHist->GetNbinsX();
  int nbinsy = inHist->GetNbinsY();

  double minx = inHist->GetXaxis()->GetBinLowEdge(1), 
         maxx = inHist->GetXaxis()->GetBinLowEdge(nbinsx + 1);

  double miny = inHist->GetYaxis()->GetBinLowEdge(1), 
         maxy = inHist->GetYaxis()->GetBinLowEdge(nbinsy + 1);
  
  auto outHist = new TH1D("ohist", "ohist", nbinsx, minx, maxx);
  for (int i = 0; i < nbinsx + 2; i++){//Including underflow and overflow
    auto h = new TH1D("hist", "hist", nbinsy, miny, maxy);
    
    for (int j = 1; j < nbinsy + 1; j++){
      double cont = inHist->GetBinContent(i,j);
      double err = inHist->GetBinError(i, j);
     
      h->SetBinContent(j, cont);
      h->SetBinError(j, err);
    }
    double fitN = 0, fitErr = 0;

    if (h->Integral() == 0){fitN = 0;}
    else{
      auto res = fitTrentoA(h);
      fitN = res[0];
      fitErr = res[1];
    }
    
    //if (fitN < 1){fitN = 0;}//Removing this at the moment
    outHist->SetBinContent(i, fitN);
    outHist->SetBinError(i, fitErr);
  }
  return outHist;
}

//Function takes in a 2D histogram with pi0 mass distribution in the y axis.
//Fits the distribution at a selected x bin and returns the integral 
//of the fit and the error on the integral. The integral range can be 
//changed by alterin the Nsig variable, which is the number of sigmas
//around the extracted pi0 mass peak to use as bouds of the integral.
vector<double>extractPi0SingleBin(TH2D* inHist, int xbin, double Nsig = n_sigma_pi0mass){
  int fint = 0;
  int nbinsy = inHist->GetNbinsY();

  double miny = inHist->GetYaxis()->GetBinLowEdge(1), 
         maxy = inHist->GetYaxis()->GetBinLowEdge(nbinsy + 1);
  
  auto h = new TH1D("hist", "hist", nbinsy, miny, maxy);
    
  for (int j = 1; j < nbinsy + 1; j++){
    double cont = inHist->GetBinContent(xbin,j);
    double err = inHist->GetBinError(xbin, j);
     
    h->SetBinContent(j, cont);
    h->SetBinError(j, err);
  }
  double fitN, fitErr;
  if (h->Integral(h->FindBin(0.1), h->FindBin(0.168)) == 0){
    fitN = 0; 
    fitErr = 0;
  }
  else{
    // trying to fit the signal starting with gaus + pol4
    // if fails it decreses pol N
    //auto fvec = findBestBack(h);
    auto fvec = findBestBack_Valerii(h);

    int fitInt = badFit(fvec,h);
    if (fitInt == 0){
      double mean  = get<2>(fvec)->GetParameter(1);
      double sigma = get<2>(fvec)->GetParameter(2);

      vector<double> rn = {mean - Nsig*sigma, mean + Nsig*sigma};
      if (fint==1) {rn = {0.1, 0.168};}
        
      fitN = get<0>(fvec)->Integral(rn[0], rn[1]) / h->GetBinWidth(1);

      auto cov = get<3>(fvec)->GetCovarianceMatrix();
      cov.ResizeTo(3,3);
      fitErr = get<1>(fvec)->IntegralError(rn[0], rn[1], get<1>(fvec)->GetParameters(), cov.GetMatrixArray());

      fitErr /= h->GetBinWidth(1);
    }
    else{
      fitN = 0;
      fitErr = 0;
    }
  }

  return (vector<double>){fitN, fitErr};
}

//Same function as above except the degree of the polynomial for the background can also
//be altered by changing polInt.
vector<double>fixextractPi0SingleBin(TH2D* inHist, int xbin, int polInt, double Nsig = n_sigma_pi0mass){
  int fint = 0;
  int nbinsy = inHist->GetNbinsY();

  double miny = inHist->GetYaxis()->GetBinLowEdge(1), 
         maxy = inHist->GetYaxis()->GetBinLowEdge(nbinsy + 1);
  
  auto h = new TH1D("hist", "hist", nbinsy, miny, maxy);
    
  for (int j = 1; j < nbinsy + 1; j++){
    double cont = inHist->GetBinContent(xbin,j);
    double err = inHist->GetBinError(xbin, j);
     
    h->SetBinContent(j, cont);
    h->SetBinError(j, err);
  }
  double fitN, fitErr;

  auto fvec = fitFunc(h, polInt);
  double mean  = get<2>(fvec)->GetParameter(1);
  double sigma = get<2>(fvec)->GetParameter(2);

  vector<double> rn = {mean - Nsig*sigma, mean + Nsig*sigma};
  if (fint==1) {rn = {0.1, 0.168};}
        
  fitN = get<0>(fvec)->Integral(rn[0], rn[1]) / h->GetBinWidth(1);

  auto cov = get<3>(fvec)->GetCovarianceMatrix();
  cov.ResizeTo(3,3);
  fitErr = get<1>(fvec)->IntegralError(rn[0], rn[1], get<1>(fvec)->GetParameters(), cov.GetMatrixArray());

  fitErr /= h->GetBinWidth(1);
  
  return (vector<double>){fitN, fitErr};
}


void PlotFitHists_sets(vector<TH1D*> inHist, vector<TF1*> signal_z, 
                                    vector<TF1*> backgr_z, vector<TF1*> sigBgr_z, const string& save_name){

  const int xPlots = 6,yPlots = 7;
  
  for (int i_set = 0; i_set < inHist.size() ; i_set+= xPlots * yPlots){
      auto c_out = new TCanvas("", "", 1800, 1800);
      c_out->Divide(xPlots, yPlots);
      for (int i_fit = i_set; i_fit < inHist.size() && i_fit < i_set + xPlots * yPlots; i_fit++){
          c_out->cd(i_fit + 1 - i_set);
          inHist.at(i_fit)->Draw("hist");
          signal_z.at(i_fit)->Draw("same");
          backgr_z.at(i_fit)->Draw("same");
          sigBgr_z.at(i_fit)->Draw("same");
      }
  
      string name_to_save = "plots/pi0fits/" + save_name + "_subsetN_" + to_string(i_set / (xPlots * yPlots))
        + ".png";
      cout << "REC fit exmp was plotted:" << name_to_save << endl;
      c_out->SaveAs(name_to_save.c_str());
  }
}


//Function takes in a 2D histogram with a pi0 mass distribution in the y axis.
//Fits the mass distribution in each x bin and returns a 1D histogram with the 
//original x axis and the bin content and error being the integral with bounds 
//of (mean - Nsig*sigma, mean + Nsig*sigma) and the error on the integral, respectively.
//Mean and sigma are taken from the gaussian portion of the fit and Nsigma can be changed.
TH1D* extractPi0(TH2D* inHist, double Nsig = n_sigma_pi0mass, bool toPlotFits = false){
  int fint = 0;
  int nbinsx = inHist->GetNbinsX();
  int nbinsy = inHist->GetNbinsY();

  double minx = inHist->GetXaxis()->GetBinLowEdge(1), 
         maxx = inHist->GetXaxis()->GetBinLowEdge(nbinsx + 1);

  double miny = inHist->GetYaxis()->GetBinLowEdge(1), 
         maxy = inHist->GetYaxis()->GetBinLowEdge(nbinsy + 1);
  
  auto outHist = new TH1D("ohist", "ohist", nbinsx, minx, maxx);

  // Saving for plotting:
  vector<TH1D*> pi0mass_zpt;
  vector<TF1*> signal_zpt;
  vector<TF1*> backgr_zpt;
  vector<TF1*> sigBgr_zpt;
  
  for (int i = 0; i < nbinsx + 2; i++){//Including underflow and overflow
    auto h = new TH1D("hist", "hist", nbinsy, miny, maxy);
      
    for (int j = 1; j < nbinsy + 1; j++){
      double cont = inHist->GetBinContent(i,j);
      double err = inHist->GetBinError(i, j);
     
      h->SetBinContent(j, cont);
      h->SetBinError(j, err);
    }

    
    double fitN, fitErr;
    // Valerii: it was equal "==0"
    if (h->Integral(h->FindBin(0.1), h->FindBin(0.168)) < 10 ){
      fitN = 0; 
      fitErr = 0;
    }
    else{
      auto fvec = findBestBack_Valerii(h);
      //auto fvec = findBestBack(h);
      int fitInt = badFit(fvec,h);
      if (fitInt == 0 /*|| fitInt == 10*/){//Changed ==0 to == 0 || ==10
        double mean  = get<2>(fvec)->GetParameter(1);
        double sigma = get<2>(fvec)->GetParameter(2);

        vector<double> rn = {mean - Nsig*sigma, mean + Nsig*sigma};
        if (fint==1) {rn = {0.1, 0.168};}
        
        fitN = get<0>(fvec)->Integral(rn[0], rn[1]) / h->GetBinWidth(1);

        auto cov = get<3>(fvec)->GetCovarianceMatrix();
        cov.ResizeTo(3,3);
        fitErr = get<1>(fvec)->IntegralError(rn[0], rn[1], get<1>(fvec)->GetParameters(), cov.GetMatrixArray());

        fitErr /= h->GetBinWidth(1);

        // drawing zpt bin fit:
        if (toPlotFits){
          int nBins_integrated = (N_Zbins) * (N_pTbins + 1) +1;
          string binName = " PT_" + to_string(i/nBins_integrated) + "_Z_" + to_string(i % nBins_integrated) ;
          auto hist_name_tmp = (string(inHist->GetName()) + binName);
          pi0mass_zpt.push_back((TH1D*)h->Clone(hist_name_tmp.c_str()));
          pi0mass_zpt.back()->SetTitle((binName).c_str());

          signal_zpt.push_back((TF1*)get<0>(fvec)->Clone());
          backgr_zpt.push_back((TF1*)get<1>(fvec)->Clone());
          sigBgr_zpt.push_back((TF1*)get<2>(fvec)->Clone());

        }
      }
      else{
        fitN = 0;
        fitErr = 0;
      }
    }

    outHist->SetBinContent(i, fitN);
    outHist->SetBinError(i, fitErr);
  }

  // Plotting is at the end because the fit "overwrote" canvas
  if (toPlotFits && pi0mass_zpt.size() > 1) {
  /*
    auto c_out = new TCanvas("", "", 2500, 2500);
    const int n_fits_hists = int(sqrt(pi0mass_zpt.size()) + 1);
    cout << "divided on:"<<n_fits_hists*(n_fits_hists-1)<<" total:"<< pi0mass_zpt.size() <<endl;
    c_out->Divide(n_fits_hists, n_fits_hists - 1);

    
    for (int i_fit = 0; i_fit < pi0mass_zpt.size(); i_fit++){
        c_out->cd(i_fit + 1);
        pi0mass_zpt.at(i_fit)->Draw("hist");
        signal_zpt.at(i_fit)->Draw("same");
        backgr_zpt.at(i_fit)->Draw("same");
        sigBgr_zpt.at(i_fit)->Draw("same");
    }

    string name_to_save = "plots/pi0fits/" + string(inHist->GetName()) + ".png";
    cout << "fit exmp was plotted:" << name_to_save << endl;
    c_out->SaveAs(name_to_save.c_str());
  */

    PlotFitHists_sets(pi0mass_zpt, signal_zpt, backgr_zpt, sigBgr_zpt, string(inHist->GetName()) );

  }
  
  return outHist;
}

void PlotFitHist(fTuple fit_funcs, TH1D* inHist,string save_name){
  auto c_out = new TCanvas("", "", 600, 600);
  inHist->Draw("hist");
  get<0>(fit_funcs)->Draw("same");
  get<1>(fit_funcs)->Draw("same");
  get<2>(fit_funcs)->Draw("same");
  string name_to_save = "plots/pi0fits/rec/" + save_name + ".png";
  c_out->SaveAs(name_to_save.c_str());
  cout << " rec fit saved:" << name_to_save << endl;
  return;
}




//This function is a higher degree version of the above function. 
//Instead of taking in a 2D histogram and returning a 1D histogram, it takes in a 3D
//histogram and returns a 2D histogram. The z axis contains the mass distributions.
TH2D* extractPi02D(TH3D* inHist, double Nsig = n_sigma_pi0mass, bool toPlotFits = false){
  int fint =0;//return to nsig after test
  int nbinsx = inHist->GetNbinsX();//pt
  int nbinsy = inHist->GetNbinsY();//z
  int nbinsz = inHist->GetNbinsZ();//pi0

  double minx = inHist->GetXaxis()->GetBinLowEdge(1), 
         maxx = inHist->GetXaxis()->GetBinLowEdge(nbinsx + 1);

  double miny = inHist->GetYaxis()->GetBinLowEdge(1), 
         maxy = inHist->GetYaxis()->GetBinLowEdge(nbinsy + 1);

  double minz = inHist->GetZaxis()->GetBinLowEdge(1), 
         maxz = inHist->GetZaxis()->GetBinLowEdge(nbinsz + 1);

    // One plot per xq2 bin
    // Storing the histograms to have access outside of the loop
    vector<TH1D*> massPi0_tmp;
    massPi0_tmp.clear();
    vector<TF1*> signal_z;
    vector<TF1*> backgr_z;
    vector<TF1*> sigBgr_z;

  
  auto outHist = new TH2D("ohist", "ohist", nbinsx, minx, maxx, nbinsy, miny, maxy);
  for (int i = 0; i < nbinsx + 2; i++){
    //cout << "pt bin:"<<i<<endl;


    
    for (int j = 0; j < nbinsy + 2; j++){
      auto h = new TH1D("hist", "hist", nbinsz, minz, maxz);
    
      for (int k = 1; k < nbinsz + 1; k++){
        double cont = inHist->GetBinContent(i, j, k);
        double err = inHist->GetBinError(i, j, k);
     
        h->SetBinContent(k, cont);
        h->SetBinError(k, err);
      }




      double fitN, fitErr;
      if (h->Integral(h->FindBin(0.1), h->FindBin(0.168)) == 0){
        fitN = 0; 
        fitErr = 0;
      }
      else{
        auto fvec = findBestBack_Valerii(h);
        //auto fvec = findBestBack(h);
        int fitInt = badFit(fvec,h);
        if (fitInt == 0/* || fitInt == 10*/){//Changed ==0 to == 0 || ==10
          double mean  = get<2>(fvec)->GetParameter(1);
          double sigma = get<2>(fvec)->GetParameter(2);

          vector<double> rn = {mean - Nsig*sigma, mean + Nsig*sigma};
          if (fint==1) {rn = {0.1, 0.168};}
        
          fitN = get<0>(fvec)->Integral(rn[0], rn[1]) / h->GetBinWidth(1);

          if (toPlotFits){
            
            auto hist_name_tmp = (string(inHist->GetName()) + " binZ:" + to_string(j));
            massPi0_tmp.push_back((TH1D*)h->Clone(hist_name_tmp.c_str()));
            massPi0_tmp.back()->SetTitle(("binZ:" + to_string(j)).c_str());
  
            signal_z.push_back((TF1*)get<0>(fvec)->Clone());
            backgr_z.push_back((TF1*)get<1>(fvec)->Clone());
            sigBgr_z.push_back((TF1*)get<2>(fvec)->Clone());

            //string plotName = string(inHist->GetName()) + "_Pt_" + to_string(i) + "_z_"+ to_string(j) ;
            //PlotFitHist(fvec, h, plotName);
          }
          //cout << "z bin:"<<j<<" fitN:"<<fitN<<endl;
          // fail here
          auto cov = get<3>(fvec)->GetCovarianceMatrix();
          cov.ResizeTo(3,3);
          fitErr = get<1>(fvec)->IntegralError(rn[0], rn[1], get<1>(fvec)->GetParameters(), cov.GetMatrixArray());

          fitErr /= h->GetBinWidth(1);
        }
        else{
          fitN = 0;
          fitErr = 0;
        }
      }

      outHist->SetBinContent(i, j, fitN);
      outHist->SetBinError(i, j, fitErr);
    }    
  }

    // PLOT EVERY xq2 BIN:
    if (toPlotFits && massPi0_tmp.size() > 1) {
      PlotFitHists_sets(massPi0_tmp, signal_z, backgr_z, sigBgr_z, string(inHist->GetName()) );
    }
  
  return outHist;
}

//Function takes in a 3D histogram with the z axis being the mass distribution
//and returns a 1D histogram of that mass distribution for a selected x and y bin.
TH1D* getMassDist(TH3D* inHist, int xbin, int ybin){
  int nbinsz  = inHist->GetNbinsZ();
  double minz = inHist->GetZaxis()->GetBinLowEdge(1), 
         maxz = inHist->GetZaxis()->GetBinLowEdge(nbinsz + 1);

  auto h = new TH1D("", "", nbinsz, minz, maxz);
  for (int i = 0; i < nbinsz + 2; i++){
    h->SetBinContent(i, inHist->GetBinContent(xbin, ybin, i));
    h->SetBinError(i, inHist->GetBinError(xbin, ybin, i));
  }
  return h;
}