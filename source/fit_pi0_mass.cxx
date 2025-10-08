// ---------------- constants-
constexpr double kPi0MeanGuess   = 0.133;   // your starting mean
constexpr double kSigmaGuess     = 0.011;   // your starting sigma

// 0.023 is to small to capture even MC
constexpr double kPeakHalfWidth  = 0.035;    // ± window for peak parameter limits
constexpr double kSidebandMin    = kPi0MeanGuess - kPeakHalfWidth;    // sideband: x < min OR x > max
constexpr double kSidebandMax    = kPi0MeanGuess + kPeakHalfWidth;
constexpr double kHardMaxFitEdge = 0.35;

constexpr double kMin_PEAKrng = 0.07;
constexpr double kMax_PEAKrng = 0.18;
constexpr double kNsig = 5.5;

// choose a reasonable upper fit edge by scanning low-yield tail
double GetLastNonTrashBinCenter(TH1D* h) {
    const double start_mass = 0.13;
    int ib = h->FindBin(start_mass);
    double last = h->GetBinContent(ib);
    // stay within safe range (avoid overflow bin)
    for (int i = ib; i <= h->GetNbinsX() - 2; ++i) {
        const double cur = h->GetBinContent(i);
        const double next = h->GetBinContent(i + 1);
        if (last + cur + next < 6) return h->GetBinCenter(i + 1);
        last = cur;
    }
    return h->GetBinCenter(h->GetNbinsX() - 1);
}

bool SkipHist_lowStat(TH1D* h){
  
  double intg_pi0mass_peak = h->Integral(h->FindBin(kPi0MeanGuess - kPeakHalfWidth), h->FindBin(kPi0MeanGuess + kPeakHalfWidth));
  double intg_next_to_pi0mass_peak = h->Integral(h->FindBin(kPi0MeanGuess + kPeakHalfWidth), h->FindBin(kPi0MeanGuess + 2*kPeakHalfWidth));

  if (h->Integral(h->FindBin(kMin_PEAKrng), h->FindBin(kMax_PEAKrng)) < 20
    || h->Integral(h->FindBin(0.1), h->FindBin(0.15)) < 20 
    || (intg_pi0mass_peak < 0.2 * intg_next_to_pi0mass_peak)
    ) return true;
  return false;
}

// fit a single Z-slice, return # pions and Error
pair<double,double> GetN_pions(TH1D* h) {
  // check how significant bckgr and if it is possible to separate them:
  double full_integral = h->Integral();
  double peak_integral = h->Integral(h->FindBin(kMin_PEAKrng), h->FindBin(kMax_PEAKrng));
  double peak_to_bckgr = (full_integral - peak_integral) / full_integral;
  
  if (peak_to_bckgr < 0.15){
    // fit with gaus only


    const double min_fit = 0.08;
    const double max_fit = 0.18;

    auto fgaus = new TF1((std::string("gausOnly_") + h->GetName()).c_str(),
                         "gaus", min_fit, max_fit);
    fgaus->SetParameters(h->GetBinContent(h->FindBin(kPi0MeanGuess)),
                         kPi0MeanGuess, kSigmaGuess);
    fgaus->SetParLimits(1, kPi0MeanGuess - 0.5 * kPeakHalfWidth,
                           kPi0MeanGuess + 0.5 * kPeakHalfWidth);
    fgaus->SetLineColor(kGreen+2);
    fgaus->SetLineWidth(2);

    auto r = h->Fit(fgaus, "QRS0");                 // quiet, use range, no draw
    h->GetListOfFunctions()->Add(fgaus);   // keep with the histogram

    // ---- draw & save PNG (includes attached function) ----
    std::filesystem::create_directories("png_only_out");
    gStyle->SetOptFit(111);

    TCanvas c((std::string("c_") + h->GetName()).c_str(), "", 900, 700);
    c.SetGrid();
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->Draw("E1");

    TLegend leg(0.60, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(h, h->GetName(), "lep");
    leg.AddEntry(fgaus, "Gaussian only", "l");
    leg.Draw();

    const std::string png = std::string("png_only_out/") + h->GetName() + "_gausOnly.png";

    double mean  = fgaus->GetParameter(1);
    double sigma = fgaus->GetParameter(2);
    vector<double> rn = {mean - kNsig*sigma, mean + kNsig*sigma};
        
    double fitN = fgaus->Integral(rn[0], rn[1]) / h->GetBinWidth(1);
    auto  cov = r->GetCovarianceMatrix();        // TMatrixDSym
    double dI = fgaus->IntegralError(rn[0], rn[1], r->GetParams(), cov.GetMatrixArray());
    double fitErr = dI / h->GetBinWidth(1);
    
    return {std::abs(fitN), std::abs(fitErr)};        // amp
  }
  else{
  
    const double edge = GetLastNonTrashBinCenter(h);
    const double max_fit = std::min(edge, kHardMaxFitEdge);

    // 1) peak prefit (gaus)
    auto fcry = new TF1((std::string("fcry_") + h->GetName()).c_str(),
                        "gaus", 0.07, std::min(0.18, edge));
    fcry->SetParameters(h->GetBinContent(h->FindBin(kPi0MeanGuess)), kPi0MeanGuess, kSigmaGuess);
    fcry->SetParLimits(1, kPi0MeanGuess - 0.5 * kPeakHalfWidth, kPi0MeanGuess + 0.5 * kPeakHalfWidth);
    fcry->SetLineColor(kGreen+2);
    fcry->SetLineWidth(2);
    h->Fit(fcry, "QRS0");

    // 2) background prefit (pol3) on sidebands
    TGraphErrors g;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double x = h->GetBinCenter(i);
        const double y = h->GetBinContent(i);
        const double ey = std::sqrt(std::max(0.0, y));
        if (x < kSidebandMin || x > kSidebandMax) {
            const int n = g.GetN();
            g.SetPoint(n, x, y);
            g.SetPointError(n, 0.0, ey);
        }
    }
    auto bkg = new TF1((std::string("bkg_") + h->GetName()).c_str(), "pol3", 0.06, max_fit);
    bkg->SetLineColor(kOrange+7);
    bkg->SetLineStyle(2);
    bkg->SetLineWidth(2);
    g.Fit(bkg, "QRS0");

    // 3) combined fit gaus(0)+pol3(3)
    auto fc = new TF1((std::string("fc_") + h->GetName()).c_str(), "gaus(0)+pol3(3)", 0.06, max_fit);
    double par[8] = {0};
    fcry->GetParameters(&par[0]);
    bkg->GetParameters(&par[3]);
    fc->SetParameters(par);
    fc->SetParLimits(1, kPi0MeanGuess - kPeakHalfWidth, kPi0MeanGuess + kPeakHalfWidth);
    fc->SetNpx(500);
    fc->SetLineColor(kRed);
    fc->SetLineWidth(2);
    auto r = h->Fit(fc, "QRS0");

    const size_t nParams_full_fit_fun = 3 + 4;
    double pf[nParams_full_fit_fun];
    for (int j = 0; j < nParams_full_fit_fun; j++){pf[j] = fc->GetParameter(j);}
    fcry->SetParameters(&pf[0]);
    bkg->SetParameters(&pf[3]);
    

    // keep functions so they’re saved with histograms (if written)
    h->GetListOfFunctions()->Add(fcry);
    h->GetListOfFunctions()->Add(bkg);
    h->GetListOfFunctions()->Add(fc);

    // Calcualte Error:
    double mean  = fc->GetParameter(1);
    double sigma = fc->GetParameter(2);
    vector<double> rn = {mean - kNsig*sigma, mean + kNsig*sigma};

    double fitN = fcry->Integral(rn[0], rn[1]) / h->GetBinWidth(1);
    auto cov = r->GetCovarianceMatrix();
    cov.ResizeTo(3,3);
    double fitErr = bkg->IntegralError(rn[0], rn[1], bkg->GetParameters(), cov.GetMatrixArray());
    fitErr /= h->GetBinWidth(1);

    // fit is unreliable so should be filtered > 0 on the next step.
    double bckgr_peak_int = bkg->Integral(rn[0], rn[1]);
    double gaus_peak_int = fcry->Integral(rn[0], rn[1]);
    if (bckgr_peak_int != 0 ){
      if ( gaus_peak_int / bckgr_peak_int < 0.35 || sigma < 0.004) return {-1,-1};
    }
        
    return {std::abs(fitN), std::abs(fitErr)};        // amp
  }
}