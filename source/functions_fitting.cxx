//using namespace clas12;
#include <string>

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
fTuple fitFunc_Valerii(TH1D *inhist, int pint = 2){

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
    const size_t nParams_full_fit_fun = 3 + (pint+1);
    double pf[nParams_full_fit_fun];
    for (int j = 0; j < nParams_full_fit_fun; j++){pf[j] = fc_poly->GetParameter(j);}
    fcry->SetParameters(&pf[0]);
    fpoly->SetParameters(&pf[3]);

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

//Function returns the z binning structure ofr the current 13 xQ2 bins
//Note that there are under flow (z < 0.2) and over flow(z>0.8) bins. This
//is for possible bin migration.
vector<double> zbins(int xqbin){
  vector<double> res;
  if (xqbin==0) {res = {0,0.2,0.275,0.350,0.425,0.500,0.575,0.650,0.725,0.8,1};}
  if (xqbin==1) {res = {0,0.2,0.227,0.257,0.281,0.317,0.359,0.419,0.515,0.8,1};}
  if (xqbin==2) {res = {0,0.2,0.239,0.275,0.311,0.353,0.401,0.461,0.557,0.8,1};}
  if (xqbin==3) {res = {0,0.2,0.257,0.299,0.341,0.389,0.437,0.497,0.581,0.8,1};}
  if (xqbin==4) {res = {0,0.2,0.287,0.335,0.377,0.419,0.467,0.521,0.593,0.8,1};}
  if (xqbin==5) {res = {0,0.2,0.335,0.377,0.413,0.449,0.485,0.527,0.581,0.8,1};}
  if (xqbin==6) {res = {0,0.2,0.227,0.257,0.293,0.329,0.377,0.437,0.539,0.8,1};}
  if (xqbin==7) {res = {0,0.2,0.239,0.275,0.317,0.359,0.407,0.473,0.563,0.8,1};}
  if (xqbin==8) {res = {0,0.2,0.269,0.311,0.353,0.401,0.449,0.503,0.581,0.8,1};}
  if (xqbin==9) {res = {0,0.2,0.305,0.353,0.389,0.425,0.467,0.509,0.569,0.8,1};}
  if (xqbin==10){res = {0,0.2,0.233,0.263,0.299,0.341,0.395,0.461,0.557,0.8,1};}
  if (xqbin==11){res = {0,0.2,0.245,0.287,0.329,0.377,0.425,0.485,0.569,0.8,1};}
  if (xqbin==12){res = {0,0.2,0.275,0.317,0.359,0.395,0.443,0.485,0.551,0.8,1};}
  if (xqbin==13){res = {0,0.2,0.245,0.281,0.317,0.359,0.407,0.461,0.539,0.8,1};}
  return res;
}




//Function returns the TH2Poly for a specific xq2-z bin.
//The structure is generated from teh z bins and the pt2 bins.
//Each z-pt2 bin is moreor less unique.
TH2Poly* makeTH2PolyZpt2(int xq2, vector<double> z){
   vector<vector<double>> p(9);
   p[8] = {0,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1};//z overflow

   if (xq2 == 0){
     p[0] = {0,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1};
     p[1] = p[0];
     p[2] = p[0];
     p[3] = p[0];
     p[4] = p[0];
     p[5] = p[0];
     p[6] = p[0];
     p[7] = p[0];
   }
   if (xq2 == 1){
     p[0] = {0,0.045,0.075,0.105,0.145,0.185,0.235,0.305,0.525};
     p[1] = {0,0.045,0.085,0.125,0.165,0.215,0.275,0.365,0.685};
     p[2] = {0,0.055,0.095,0.135,0.185,0.245,0.315,0.425,0.805};
     p[3] = {0,0.055,0.105,0.155,0.205,0.275,0.355,0.475,1};
     p[4] = {0,0.055,0.115,0.165,0.215,0.285,0.385,0.525,1};
     p[5] = {0,0.065,0.125,0.185,0.245,0.315,0.415,0.565,1};
     p[6] = {0,0.075,0.145,0.205,0.265,0.345,0.445,0.605,1};
     p[7] = {0,0.075,0.145,0.215,0.285,0.375,0.485,0.645,1};
   }
   if (xq2==2){
     p[0] = {0,0.035,0.055,0.085,0.115,0.155,0.195,0.255,0.525};
     p[1] = {0,0.035,0.065,0.105,0.145,0.185,0.245,0.325,0.725};
     p[2] = {0,0.035,0.075,0.115,0.165,0.215,0.285,0.385,0.845};
     p[3] = {0,0.045,0.085,0.135,0.185,0.235,0.315,0.435,1};
     p[4] = {0,0.045,0.095,0.145,0.195,0.265,0.345,0.485,1};
     p[5] = {0,0.045,0.095,0.155,0.215,0.285,0.375,0.515,1};
     p[6] = {0,0.045,0.105,0.165,0.225,0.305,0.405,0.555,1};
     p[7] = {0,0.035,0.085,0.145,0.225,0.305,0.405,0.565,1};
   }
   if (xq2==3){
     p[0] = {0,0.025,0.045,0.065,0.095,0.125,0.165,0.225,0.475};
     p[1] = {0,0.025,0.055,0.095,0.125,0.165,0.215,0.285,0.656};
     p[2] = {0,0.025,0.065,0.105,0.145,0.195,0.255,0.345,0.785};
     p[3] = {0,0.035,0.065,0.115,0.165,0.215,0.285,0.395,1};
     p[4] = {0,0.035,0.075,0.115,0.175,0.235,0.305,0.435,1};
     p[5] = {0,0.035,0.075,0.125,0.175,0.245,0.335,0.465,1};
     p[6] = {0,0.025,0.065,0.125,0.185,0.265,0.355,0.495,1};
     p[7] = {0,0.025,0.055,0.095,0.155,0.225,0.325,0.455,1};
   } 
   if (xq2==4){
     p[0] = {0,0.015,0.035,0.055,0.085,0.115,0.155,0.205,0.485};
     p[1] = {0,0.025,0.045,0.075,0.105,0.145,0.195,0.275,0.675};
     p[2] = {0,0.025,0.055,0.085,0.115,0.165,0.225,0.315,0.855};
     p[3] = {0,0.025,0.055,0.085,0.135,0.185,0.255,0.355,1};
     p[4] = {0,0.025,0.055,0.095,0.135,0.195,0.265,0.385,1};
     p[5] = {0,0.025,0.055,0.095,0.145,0.205,0.275,0.395,1};
     p[6] = {0,0.025,0.045,0.085,0.135,0.195,0.285,0.395,1};
     p[7] = {0,0.015,0.035,0.065,0.095,0.155,0.225,0.335,1};
   }
   if (xq2==5){
     p[0] = {0,0.015,0.025,0.045,0.065,0.105,0.145,0.205,0.535};
     p[1] = {0,0.015,0.035,0.055,0.075,0.115,0.165,0.235,0.625};
     p[2] = {0,0.015,0.035,0.055,0.075,0.115,0.165,0.245,0.835};
     p[3] = {0,0.015,0.035,0.055,0.085,0.115,0.165,0.255,0.805};
     p[4] = {0,0.015,0.025,0.055,0.075,0.115,0.165,0.255,0.815};
     p[5] = {0,0.015,0.025,0.045,0.075,0.115,0.165,0.255,0.815};
     p[6] = {0,0.015,0.025,0.045,0.075,0.105,0.155,0.245,0.715};
     p[7] = {0,0.015,0.025,0.035,0.055,0.085,0.125,0.205,0.695};
   }
   if (xq2==6){
     p[0] = {0,0.035,0.065,0.105,0.135,0.175,0.215,0.285,0.495};
     p[1] = {0,0.045,0.085,0.115,0.155,0.205,0.265,0.345,0.615};
     p[2] = {0,0.045,0.085,0.135,0.175,0.235,0.305,0.405,0.805};
     p[3] = {0,0.045,0.095,0.145,0.195,0.255,0.345,0.475,1};
     p[4] = {0,0.055,0.105,0.155,0.215,0.285,0.375,0.515,1};
     p[5] = {0,0.055,0.115,0.165,0.235,0.305,0.395,0.555,1};
     p[6] = {0,0.065,0.115,0.175,0.245,0.325,0.425,0.585,1};
     p[7] = {0,0.055,0.115,0.175,0.255,0.335,0.445,0.605,1};
   }
   if(xq2==7){
     p[0] = {0,0.025,0.055,0.085,0.115,0.145,0.185,0.245,0.495};
     p[1] = {0,0.035,0.065,0.095,0.135,0.175,0.225,0.305,0.615};
     p[2] = {0,0.035,0.065,0.115,0.155,0.205,0.265,0.375,0.865};
     p[3] = {0,0.035,0.075,0.125,0.175,0.225,0.305,0.425,1};
     p[4] = {0,0.035,0.075,0.125,0.185,0.245,0.335,0.465,1};
     p[5] = {0,0.035,0.085,0.135,0.195,0.265,0.355,0.495,1};
     p[6] = {0,0.035,0.085,0.135,0.205,0.275,0.375,0.525,1};
     p[7] = {0,0.025,0.065,0.115,0.185,0.265,0.365,0.515,1};
   }
   if(xq2==8){
     p[0] = {0,0.015,0.035,0.065,0.085,0.115,0.155,0.215,0.495};
     p[1] = {0,0.025,0.045,0.075,0.115,0.155,0.205,0.275,0.635};
     p[2] = {0,0.025,0.055,0.085,0.125,0.175,0.235,0.325,0.805};
     p[3] = {0,0.025,0.055,0.085,0.135,0.185,0.255,0.365,1};
     p[4] = {0,0.025,0.055,0.095,0.145,0.205,0.275,0.385,1};
     p[5] = {0,0.025,0.055,0.095,0.145,0.205,0.285,0.405,1};
     p[6] = {0,0.025,0.055,0.095,0.145,0.205,0.295,0.425,1};
     p[7] = {0,0.015,0.045,0.075,0.115,0.165,0.245,0.365,1};
   }
   if(xq2==9){
     p[0] = {0,0.015,0.025,0.045,0.065,0.095,0.135,0.185,0.435};
     p[1] = {0,0.015,0.035,0.055,0.075,0.105,0.155,0.225,0.615};
     p[2] = {0,0.015,0.035,0.055,0.085,0.115,0.165,0.245,0.735};
     p[3] = {0,0.015,0.035,0.055,0.085,0.115,0.175,0.255,0.835};
     p[4] = {0,0.015,0.035,0.055,0.085,0.125,0.175,0.265,1};
     p[5] = {0,0.015,0.035,0.055,0.085,0.115,0.165,0.255,1};
     p[6] = {0,0.015,0.035,0.045,0.075,0.115,0.165,0.255,0.855};
     p[7] = {0,0.015,0.025,0.045,0.065,0.085,0.135,0.205,0.695};
   }
   if(xq2==10){
     p[0] = {0,0.025,0.055,0.085,0.125,0.155,0.205,0.265,0.485};
     p[1] = {0,0.035,0.065,0.105,0.145,0.185,0.245,0.325,0.605};
     p[2] = {0,0.035,0.075,0.115,0.165,0.215,0.285,0.385,0.815};
     p[3] = {0,0.045,0.085,0.125,0.175,0.245,0.325,0.445,1};
     p[4] = {0,0.045,0.085,0.135,0.195,0.255,0.345,0.495,1};
     p[5] = {0,0.045,0.095,0.145,0.205,0.285,0.375,0.525,1};
     p[6] = {0,0.045,0.085,0.145,0.205,0.285,0.385,0.535,1};
     p[7] = {0,0.035,0.085,0.135,0.195,0.285,0.385,0.535,1};
   }
   if(xq2==11){
     p[0] = {0,0.025,0.045,0.065,0.095,0.125,0.165,0.215,0.455};
     p[1] = {0,0.025,0.045,0.075,0.115,0.155,0.205,0.275,0.585};
     p[2] = {0,0.025,0.055,0.095,0.135,0.185,0.245,0.335,0.805};
     p[3] = {0,0.025,0.065,0.095,0.145,0.205,0.275,0.385,1};
     p[4] = {0,0.025,0.065,0.105,0.155,0.215,0.295,0.415,1};
     p[5] = {0,0.025,0.065,0.105,0.155,0.225,0.305,0.435,1};
     p[6] = {0,0.025,0.055,0.095,0.155,0.215,0.315,0.445,1};
     p[7] = {0,0.025,0.045,0.085,0.125,0.185,0.265,0.395,1};
   }
   if(xq2==12){
     p[0] = {0,0.015,0.025,0.045,0.065,0.095,0.125,0.175,0.415};
     p[1] = {0,0.015,0.035,0.055,0.075,0.115,0.155,0.225,0.575};
     p[2] = {0,0.015,0.035,0.055,0.085,0.125,0.175,0.255,0.735};
     p[3] = {0,0.015,0.035,0.065,0.095,0.135,0.185,0.275,0.845};
     p[4] = {0,0.015,0.035,0.065,0.095,0.135,0.195,0.285,1};
     p[5] = {0,0.015,0.035,0.065,0.095,0.135,0.195,0.295,1};
     p[6] = {0,0.015,0.035,0.055,0.085,0.125,0.185,0.275,1};
     p[7] = {0,0.015,0.035,0.045,0.075,0.105,0.155,0.245,0.825};
   }
   if(xq2==13){
     p[0] = {0,0.015,0.035,0.065,0.085,0.115,0.155,0.205,0.445};
     p[1] = {0,0.025,0.045,0.075,0.105,0.145,0.195,0.265,0.615};
     p[2] = {0,0.025,0.045,0.075,0.115,0.165,0.225,0.315,0.755};
     p[3] = {0,0.025,0.055,0.085,0.135,0.185,0.245,0.355,1};
     p[4] = {0,0.025,0.055,0.085,0.135,0.185,0.255,0.375,1};
     p[5] = {0,0.025,0.055,0.085,0.125,0.185,0.265,0.375,1};
     p[6] = {0,0.025,0.045,0.085,0.125,0.185,0.255,0.385,1};
     p[7] = {0,0.025,0.045,0.075,0.115,0.165,0.245,0.365,1};
   }
   
   //pt overflow bin; pt > ptmx in z bin
   for (int i = 0; i < p.size(); i++){
     p[i].push_back(1.5);
   }

   auto poly = new TH2Poly(); poly->Sumw2();
   poly->SetLineStyle(1);
   poly->SetLineWidth(2);

   for (int i = 0; i < z.size()-1; i++){
     if (i == 0){
       for (int j = 0; j < p[i].size()-1; j++){
         poly->AddBin(z[i], p[8][j], z[i+1], p[8][j+1]);
       }
     }
     else{
       for (int j = 0; j < p[i-1].size()-1; j++){
         poly->AddBin(z[i], p[i-1][j], z[i+1], p[i-1][j+1]);
       }
     }
   }

    return poly;
 }




//Function returns a 2D vector representing the pt2 bin edges for a 
//chosen xQ2 bin

// Commented by Valerii to not change the name:
/*
vector<vector<double>> pbins(int &xq2){
   vector<vector<double>> p(9);
   p[8] = {0,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1};//z overflow

   auto z = zbins(xq2);
   if (xq2 == 0){
     p[0] = {0,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1};
     p[1] = p[0];
     p[2] = p[0];
     p[3] = p[0];
     p[4] = p[0];
     p[5] = p[0];
     p[6] = p[0];
     p[7] = p[0];
   }
   if (xq2 == 1){
     p[0] = {0,0.045,0.075,0.105,0.145,0.185,0.235,0.305,0.525};
     p[1] = {0,0.045,0.085,0.125,0.165,0.215,0.275,0.365,0.685};
     p[2] = {0,0.055,0.095,0.135,0.185,0.245,0.315,0.425,0.805};
     p[3] = {0,0.055,0.105,0.155,0.205,0.275,0.355,0.475,1};
     p[4] = {0,0.055,0.115,0.165,0.215,0.285,0.385,0.525,1};
     p[5] = {0,0.065,0.125,0.185,0.245,0.315,0.415,0.565,1};
     p[6] = {0,0.075,0.145,0.205,0.265,0.345,0.445,0.605,1};
     p[7] = {0,0.075,0.145,0.215,0.285,0.375,0.485,0.645,1};
   }
   if (xq2==2){
     p[0] = {0,0.035,0.055,0.085,0.115,0.155,0.195,0.255,0.525};
     p[1] = {0,0.035,0.065,0.105,0.145,0.185,0.245,0.325,0.725};
     p[2] = {0,0.035,0.075,0.115,0.165,0.215,0.285,0.385,0.845};
     p[3] = {0,0.045,0.085,0.135,0.185,0.235,0.315,0.435,1};
     p[4] = {0,0.045,0.095,0.145,0.195,0.265,0.345,0.485,1};
     p[5] = {0,0.045,0.095,0.155,0.215,0.285,0.375,0.515,1};
     p[6] = {0,0.045,0.105,0.165,0.225,0.305,0.405,0.555,1};
     p[7] = {0,0.035,0.085,0.145,0.225,0.305,0.405,0.565,1};
   }
   if (xq2==3){
     p[0] = {0,0.025,0.045,0.065,0.095,0.125,0.165,0.225,0.475};
     p[1] = {0,0.025,0.055,0.095,0.125,0.165,0.215,0.285,0.656};
     p[2] = {0,0.025,0.065,0.105,0.145,0.195,0.255,0.345,0.785};
     p[3] = {0,0.035,0.065,0.115,0.165,0.215,0.285,0.395,1};
     p[4] = {0,0.035,0.075,0.115,0.175,0.235,0.305,0.435,1};
     p[5] = {0,0.035,0.075,0.125,0.175,0.245,0.335,0.465,1};
     p[6] = {0,0.025,0.065,0.125,0.185,0.265,0.355,0.495,1};
     p[7] = {0,0.025,0.055,0.095,0.155,0.225,0.325,0.455,1};
   }
   if (xq2==4){
     p[0] = {0,0.015,0.035,0.055,0.085,0.115,0.155,0.205,0.485};
     p[1] = {0,0.025,0.045,0.075,0.105,0.145,0.195,0.275,0.675};
     p[2] = {0,0.025,0.055,0.085,0.115,0.165,0.225,0.315,0.855};
     p[3] = {0,0.025,0.055,0.085,0.135,0.185,0.255,0.355,1};
     p[4] = {0,0.025,0.055,0.095,0.135,0.195,0.265,0.385,1};
     p[5] = {0,0.025,0.055,0.095,0.145,0.205,0.275,0.395,1};
     p[6] = {0,0.025,0.045,0.085,0.135,0.195,0.285,0.395,1};
     p[7] = {0,0.015,0.035,0.065,0.095,0.155,0.225,0.335,1};
   }
   if (xq2==5){
     p[0] = {0,0.015,0.025,0.045,0.065,0.105,0.145,0.205,0.535};
     p[1] = {0,0.015,0.035,0.055,0.075,0.115,0.165,0.235,0.625};
     p[2] = {0,0.015,0.035,0.055,0.075,0.115,0.165,0.245,0.835};
     p[3] = {0,0.015,0.035,0.055,0.085,0.115,0.165,0.255,0.805};
     p[4] = {0,0.015,0.025,0.055,0.075,0.115,0.165,0.255,0.815};
     p[5] = {0,0.015,0.025,0.045,0.075,0.115,0.165,0.255,0.815};
     p[6] = {0,0.015,0.025,0.045,0.075,0.105,0.155,0.245,0.715};
     p[7] = {0,0.015,0.025,0.035,0.055,0.085,0.125,0.205,0.695};
   }
   if (xq2==6){
     p[0] = {0,0.035,0.065,0.105,0.135,0.175,0.215,0.285,0.495};
     p[1] = {0,0.045,0.085,0.115,0.155,0.205,0.265,0.345,0.615};
     p[2] = {0,0.045,0.085,0.135,0.175,0.235,0.305,0.405,0.805};
     p[3] = {0,0.045,0.095,0.145,0.195,0.255,0.345,0.475,1};
     p[4] = {0,0.055,0.105,0.155,0.215,0.285,0.375,0.515,1};
     p[5] = {0,0.055,0.115,0.165,0.235,0.305,0.395,0.555,1};
     p[6] = {0,0.065,0.115,0.175,0.245,0.325,0.425,0.585,1};
     p[7] = {0,0.055,0.115,0.175,0.255,0.335,0.445,0.605,1};
   }
   if(xq2==7){
     p[0] = {0,0.025,0.055,0.085,0.115,0.145,0.185,0.245,0.495};
     p[1] = {0,0.035,0.065,0.095,0.135,0.175,0.225,0.305,0.615};
     p[2] = {0,0.035,0.065,0.115,0.155,0.205,0.265,0.375,0.865};
     p[3] = {0,0.035,0.075,0.125,0.175,0.225,0.305,0.425,1};
     p[4] = {0,0.035,0.075,0.125,0.185,0.245,0.335,0.465,1};
     p[5] = {0,0.035,0.085,0.135,0.195,0.265,0.355,0.495,1};
     p[6] = {0,0.035,0.085,0.135,0.205,0.275,0.375,0.525,1};
     p[7] = {0,0.025,0.065,0.115,0.185,0.265,0.365,0.515,1};
   }
   if(xq2==8){
     p[0] = {0,0.015,0.035,0.065,0.085,0.115,0.155,0.215,0.495};
     p[1] = {0,0.025,0.045,0.075,0.115,0.155,0.205,0.275,0.635};
     p[2] = {0,0.025,0.055,0.085,0.125,0.175,0.235,0.325,0.805};
     p[3] = {0,0.025,0.055,0.085,0.135,0.185,0.255,0.365,1};
     p[4] = {0,0.025,0.055,0.095,0.145,0.205,0.275,0.385,1};
     p[5] = {0,0.025,0.055,0.095,0.145,0.205,0.285,0.405,1};
     p[6] = {0,0.025,0.055,0.095,0.145,0.205,0.295,0.425,1};
     p[7] = {0,0.015,0.045,0.075,0.115,0.165,0.245,0.365,1};
   }

   if(xq2==9){
     p[0] = {0,0.015,0.025,0.045,0.065,0.095,0.135,0.185,0.435};
     p[1] = {0,0.015,0.035,0.055,0.075,0.105,0.155,0.225,0.615};
     p[2] = {0,0.015,0.035,0.055,0.085,0.115,0.165,0.245,0.735};
     p[3] = {0,0.015,0.035,0.055,0.085,0.115,0.175,0.255,0.835};
     p[4] = {0,0.015,0.035,0.055,0.085,0.125,0.175,0.265,1};
     p[5] = {0,0.015,0.035,0.055,0.085,0.115,0.165,0.255,1};
     p[6] = {0,0.015,0.035,0.045,0.075,0.115,0.165,0.255,0.855};
     p[7] = {0,0.015,0.025,0.045,0.065,0.085,0.135,0.205,0.695};
   }
   if(xq2==10){
     p[0] = {0,0.025,0.055,0.085,0.125,0.155,0.205,0.265,0.485};
     p[1] = {0,0.035,0.065,0.105,0.145,0.185,0.245,0.325,0.605};
     p[2] = {0,0.035,0.075,0.115,0.165,0.215,0.285,0.385,0.815};
     p[3] = {0,0.045,0.085,0.125,0.175,0.245,0.325,0.445,1};
     p[4] = {0,0.045,0.085,0.135,0.195,0.255,0.345,0.495,1};
     p[5] = {0,0.045,0.095,0.145,0.205,0.285,0.375,0.525,1};
     p[6] = {0,0.045,0.085,0.145,0.205,0.285,0.385,0.535,1};
     p[7] = {0,0.035,0.085,0.135,0.195,0.285,0.385,0.535,1};
   }
   if(xq2==11){
     p[0] = {0,0.025,0.045,0.065,0.095,0.125,0.165,0.215,0.455};
     p[1] = {0,0.025,0.045,0.075,0.115,0.155,0.205,0.275,0.585};
     p[2] = {0,0.025,0.055,0.095,0.135,0.185,0.245,0.335,0.805};
     p[3] = {0,0.025,0.065,0.095,0.145,0.205,0.275,0.385,1};
     p[4] = {0,0.025,0.065,0.105,0.155,0.215,0.295,0.415,1};
     p[5] = {0,0.025,0.065,0.105,0.155,0.225,0.305,0.435,1};
     p[6] = {0,0.025,0.055,0.095,0.155,0.215,0.315,0.445,1};
     p[7] = {0,0.025,0.045,0.085,0.125,0.185,0.265,0.395,1};
   }
   if(xq2==12){
     p[0] = {0,0.015,0.025,0.045,0.065,0.095,0.125,0.175,0.415};
     p[1] = {0,0.015,0.035,0.055,0.075,0.115,0.155,0.225,0.575};
     p[2] = {0,0.015,0.035,0.055,0.085,0.125,0.175,0.255,0.735};
     p[3] = {0,0.015,0.035,0.065,0.095,0.135,0.185,0.275,0.845};
     p[4] = {0,0.015,0.035,0.065,0.095,0.135,0.195,0.285,1};
     p[5] = {0,0.015,0.035,0.065,0.095,0.135,0.195,0.295,1};
     p[6] = {0,0.015,0.035,0.055,0.085,0.125,0.185,0.275,1};
     p[7] = {0,0.015,0.035,0.045,0.075,0.105,0.155,0.245,0.825};
   }
   if(xq2==13){
     p[0] = {0,0.015,0.035,0.065,0.085,0.115,0.155,0.205,0.445};
     p[1] = {0,0.025,0.045,0.075,0.105,0.145,0.195,0.265,0.615};
     p[2] = {0,0.025,0.045,0.075,0.115,0.165,0.225,0.315,0.755};
     p[3] = {0,0.025,0.055,0.085,0.135,0.185,0.245,0.355,1};
     p[4] = {0,0.025,0.055,0.085,0.135,0.185,0.255,0.375,1};
     p[5] = {0,0.025,0.055,0.085,0.125,0.185,0.265,0.375,1};
     p[6] = {0,0.025,0.045,0.085,0.125,0.185,0.255,0.385,1};
     p[7] = {0,0.025,0.045,0.075,0.115,0.165,0.245,0.365,1};
   }
   
   //pt overflow bin; pt > ptmx in z bin
   for (int i = 0; i < p.size(); i++){
     p[i].push_back(1.5);
   }

   vector<vector<double>> pf(p.size() + 1);
   for(int i = 0; i < p.size(); i++){
     pf[i+1] = p[i]; 
   }
   pf[0] = p[8];
   return pf;
}
*/


//Function returns the 2D z-pt2 bin width of a selected xQ2-z-pt2 bin
double zpt2BinWidth(int xq2, int inZ, int inP){
  //auto z = zbins(xq2);
  auto z = zbins_Valerii(xq2);
  auto p = pbins(xq2);

  return (z[inZ+1] - z[inZ]) * (p[inZ][inP+1] - p[inZ][inP]);
}

//Function returns the z bin with for a selected xQ2-z bin
double zBinWidth(int xq2, int inZ){
  //auto z = zbins(xq2);
  auto z = zbins_Valerii(xq2);
  return z[inZ+1] - z[inZ];
}

//Function returns the pt2 bin with for a selected xQ2-z-pt2 bin
double pt2BinWidth(int xq2, int inZ, int inP){
  auto p = pbins(xq2);

  return p[inZ][inP+1] - p[inZ][inP];
}

//Function returns the bin center of a z bin using the left and right bin edges
double zBinCenter(int xq2, int inZ){
  //auto z = zbins(xq2);
  auto z = zbins_Valerii(xq2);
  return (z[inZ+1] + z[inZ])*0.5;
}

//Function returns the bin center of a pt2 bin using the left and right bin edges
double pt2BinCenter(int xq2, int inZ, int inP){
  auto p = pbins(xq2);

  return (p[inZ][inP+1] + p[inZ][inP])*0.5;
}

//Funtion returns the distance between the bincenter and the left and right
//edges for asymmetric bin errors
vector<double> zBinEdges(int xq2, int inZ, double binCenter){
  //auto z = zbins(xq2);
  auto z = zbins_Valerii(xq2);
  return (vector<double>){binCenter - z[inZ], z[inZ+1] - binCenter};
}

//Funtion returns the distance between the bincenter and the left and right
//edges for asymmetric bin errors
vector<double> pt2BinEdges(int xq2, int inZ, int inP, double binCenter){
  auto p = pbins(xq2);
  return (vector<double>){binCenter - p[inZ][inP], p[inZ][inP+1] - binCenter};
}

//Function takes in a 2D histogram and scales each bin by 1/binwidth
//the bin width is the 2D z-pt2 binwidth. This function assumes that 
//the histogram has z on the x axis and pt2 on the y axis.
void scaleByzpt2BinWidth(TH2D* inhist, int xq2){
  for (int i = 1; i < inhist->GetNbinsX() + 1; i++){
   for (int j = 1; j < inhist->GetNbinsY() + 1; j++){
     double cont = inhist->GetBinContent(i, j);
     double err  = inhist->GetBinError(i, j);
     //cout<<"my guess the problem is there:"<<endl;
     cont /= zpt2BinWidth(xq2, i-1, j-1);
     err  /= zpt2BinWidth(xq2, i-1, j-1);
     inhist->SetBinContent(i, j, cont);
     inhist->SetBinError(i, j, err);
   }
  }
}

//Function takes in a 2D histogram and scales each bin by binwidth
//the bin width is the pt2 binwidth. This function assumes that 
//the histogram has z on the x axis and pt2 on the y axis.
void scaleBypt2BinWidth(TH2D* inhist, int xq2){
  for (int i = 1; i < inhist->GetNbinsX() + 1; i++){
   for (int j = 1; j < inhist->GetNbinsY() + 1; j++){
     double cont = inhist->GetBinContent(i, j);
     double err  = inhist->GetBinError(i, j);
     cont *= pt2BinWidth(xq2, i-1, j-1);
     err  *= pt2BinWidth(xq2, i-1, j-1);
     inhist->SetBinContent(i, j, cont);
     inhist->SetBinError(i, j, err);
   }
  }
}

//Function generates a 2D histogram
  TH2Poly* makeTH2Poly(){
    auto poly = new TH2Poly();
    double x1[] = {0.0835, 0.15,  0.15, 0.12};
    double y1[] = {1.3,    2.28,  1.38,  1.3};
    double x2[] = {0.15,   0.15,  0.24, 0.24, 0.20};
    double y2[] = {1.38,   1.98,  2.75,  1.5, 1.45};
    double x3[] = {0.15,   0.15,  0.24, 0.24};
    double y3[] = {1.98,   2.28, 3.625, 2.75};
    double x4[] = {0.24,   0.24,  0.34, 0.34, 0.30, 0.27};
    double y4[] = {1.5,    2.75,  3.63,  1.6, 1.56, 1.53};
    double x5[] = {0.24,   0.24,  0.34, 0.34};
    double y5[] = {2.75,  3.625,  5.12, 3.63};
    double x6[] = {0.34,   0.34,  0.45, 0.45};
    double y6[] = {1.6,    3.63,  4.7,  2.52};
    double x7[] = {0.34,   0.34,  0.45, 0.45};
    double y7[] = {3.63,   5.12,  6.76,  4.7};
    double x8[] = {0.45,   0.45, 0.708, 0.64, 0.57, 0.50};
    double y8[] = {2.52,    4.7,  7.42,  5.4, 4.05, 3.05};
    double x9[] = {0.45,   0.45,  0.677, 0.7896, 0.75, 0.708};
    double y9[] = {4.7,    6.76, 10.185, 11.351, 9.52,  7.42};

    poly->AddBin(4, x1, y1);
    poly->AddBin(5, x2, y2);
    poly->AddBin(4, x3, y3);
    poly->AddBin(5, x4, y4);
    poly->AddBin(4, x5, y5);
    poly->AddBin(4, x6, y6);
    poly->AddBin(4, x7, y7);
    poly->AddBin(6, x8, y8);
    poly->AddBin(6, x9, y9);

    return poly;
 }


 //Function makes and 2D histigram of Richard's original binning scheme
 TH2Poly* makeTH2PolyRich(){
    auto poly = new TH2Poly();
    // X axis is Xb variable
    // Y axis is Q2 variable
    double x1[] = {0.126602, 0.15,  0.24,    0.24,   0.15,     0.126602};
    double y1[] = {2,        2.28,  3.625,   2.75,   2,        2};
    double x2[] = {0.15,     0.24,  0.24,    0.15};
    double y2[] = {2,        2.75,  2,       2};
    double x3[] = {0.24,     0.24,  0.34,    0.34,   0.24};
    double y3[] = {2.75,     3.625, 5.12,    3.63,   2.75};
    double x4[] = {0.24,     0.24,  0.34,    0.34,   0.24};
    double y4[] = {2,        2.75,  3.63,    2, 2};
    double x5[] = {0.34,     0.34,  0.45,    0.45,   0.34};
    double y5[] = {3.63,     5.12,  6.76,    4.7,    3.63};
    double x6[] = {0.34,     0.34,  0.45,    0.45,   0.387826, 0.34};
    double y6[] = {2,        3.63,  4.7,     2.52,   2, 2};
    double x7[] = {0.45,     0.45,  0.677,   0.7896, 0.75,     0.708, 0.45};
    double y7[] = {4.7,      6.76,  10.185,  11.351, 9.52,     7.42,  4.7};
    double x8[] = {0.45,     0.45,  0.708,   0.64,   0.57,     0.50,  0.45};
    double y8[] = {2.52,     4.7,   7.42,    5.4,    4.05,     3.05,  2.52};

    poly->AddBin(6, x1, y1);
    poly->AddBin(4, x2, y2);
    poly->AddBin(5, x3, y3);
    poly->AddBin(5, x4, y4);
    poly->AddBin(5, x5, y5);
    poly->AddBin(6, x6, y6);
    poly->AddBin(7, x7, y7);
    poly->AddBin(7, x8, y8);

    poly->SetLineStyle(1);
    poly->SetLineWidth(10);

    return poly;
 }

 //Fuction geenrates a TLine vector of the outline for Richard's original binning
 vector<TLine> makeRichOutline(int bin = 0){
    vector<TLine> tvec(26);
    
    //bin 1
    tvec[0] = TLine(0.126602, 2, 0.15, 2); 
    tvec[1] = TLine(0.15, 2, 0.24, 2.75); 
    tvec[2] = TLine(0.24, 2.75, 0.24, 3.625);     
    tvec[3] = TLine(0.24, 3.625, 0.15, 2.28); 
    tvec[4] = TLine(0.15, 2.28, 0.126602, 2);   

    //bin 2
    tvec[5] = TLine(0.15, 2, 0.24, 2); 
    tvec[6] = TLine(0.24, 2, 0.24, 2.75); 

    //bin 3
    tvec[7] = TLine(0.24, 3.625, 0.34, 5.12); 
    tvec[8] = TLine(0.34, 5.12, 0.34, 3.63); 
    tvec[9] = TLine(0.34, 3.63, 0.24, 2.75); 

    //bin 4
    tvec[10] = TLine(0.24, 2, 0.34, 2); 
    tvec[11] = TLine(0.34, 2, 0.34, 3.63); 

    //bin 5
    tvec[12] = TLine(0.34, 5.12, 0.45, 6.76); 
    tvec[13] = TLine(0.45, 6.76, 0.45, 4.7); 
    tvec[14] = TLine(0.45, 4.7, 0.34, 3.63); 

    //bin 6
    tvec[15] = TLine(0.34, 2, 0.387826, 2); 
    tvec[16] = TLine(0.387826, 2, 0.45, 2.52); 
    tvec[17] = TLine(0.45, 2.52, 0.45, 4.7); 
  
    //bin 7 
    tvec[18] = TLine(0.45, 6.76, 0.677, 10.185); 
    tvec[19] = TLine(0.677, 10.185, 0.7896, 11.351); 
    tvec[20] = TLine(0.7896, 11.351, 0.708, 7.42);    
    tvec[21] = TLine(0.708, 7.42, 0.45, 4.7);   

    //bin 8 
    tvec[22] = TLine(0.708, 7.42, 0.64, 5.4); 
    tvec[23] = TLine(0.64, 5.4, 0.57, 4.05); 
    tvec[24] = TLine(0.57, 4.05, 0.5, 3.05);    
    tvec[25] = TLine(0.5, 3.05, 0.45, 2.52);  

    for (int i = 0; i < tvec.size(); i++){
      tvec[i].SetLineWidth(2);
    }

    vector<vector<int>> blin(9);
    blin[1] = {0, 1, 2, 3, 4};
    blin[2] = {1, 5, 6};
    blin[3] = {7, 8, 9, 2};
    blin[4] = {10, 11, 6, 9};
    blin[5] = {12, 13, 14, 8};
    blin[6] = {15, 16, 17, 14, 11};
    blin[7] = {18, 19, 20, 21, 13};
    blin[8] = {22, 23, 24, 25, 21, 17};

    if (bin > 0 && bin < 9){
      for (int i = 0; i < blin[bin].size(); i++){
        tvec[blin[bin][i]].SetLineWidth(4);
        tvec[blin[bin][i]].SetLineColor(kRed);
      }
    }

    return tvec;
 }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


 auto makeTH2PolyMars(){
    auto poly = new TH2Poly();
//1st row
    double x1[] = {0.126,  0.18,    0.18,     0.126};//0.18
    double y1[] = {2,      2.658,   2,        2};

    double x2[] = {0.18,   0.18,    0.21,     0.21,   0.18};//.18
    double y2[] = {2,      2.658,   2.658,    2,      2};

    double x3[] = {0.21,   0.21,    0.24,     0.24,   0.21};//.18
    double y3[] = {2,      2.658,   2.658,    2,      2};

    double x4[] = {0.24,   0.28,    0.28,     0.24};
    double y4[] = {2.658,  2.658,   2.00,     2};

    double x5[] = {0.28,   0.28,    0.45,     0.45,   0.387826, 0.28};
    double y5[] = {2.00,   2.658,   2.658,    2.658,  2,        2};
//2nd row
    double x6[] = {0.18,   0.24,    0.24,     0.18};//.18
    double y6[] = {2.658,  3.625,   2.658,    2.658};

    double x7[] = {0.24,   0.24,    0.28,     0.28,   0.24};
    double y7[] = {2.658,  3.625,   3.625,    2.658,  2.658};
   
    double x8[] = {0.28,   0.28,    0.34,     0.34,   0.28};
    double y8[] = {2.658,  3.625,   3.625,    2.658,  2.658};

    double x9[] = {0.34,   0.34,    0.54025,  0.45,   0.34};
    double y9[] = {2.658,  3.625,   3.625,    2.658,  2.658};
//3rd row
    double x10[] = {0.24,  0.34,    0.34,     0.24};
    double y10[] = {3.625, 5.12,    3.625,    3.625};

    double x11[] = {0.34,  0.34,    0.40,     0.40,   0.34};
    double y11[] = {3.625, 5.12,    5.12,     3.625,  3.625};

    double x12[] = {0.40,  0.40,    0.6234,   0.57,   0.54025,  0.40};
    double y12[] = {3.625, 5.12,    5.12,     4.05,   3.625,    3.625};
//4th row
    double x13[] = {0.34,  0.677,   0.7896,   0.75,   0.708,    0.64,  0.6234, 0.34};
    double y13[] = {5.12,  10.185,  11.351,   9.52,   7.42,     5.4,   5.12,   5.12};

    poly->AddBin(4, x1, y1);
    poly->AddBin(5, x2, y2);
    poly->AddBin(5, x3, y3);
    poly->AddBin(4, x4, y4);
    poly->AddBin(6, x5, y5);
    poly->AddBin(4, x6, y6);
    poly->AddBin(5, x7, y7);
    poly->AddBin(5, x8, y8);
    poly->AddBin(5, x9, y9);
    poly->AddBin(4, x10, y10);
    poly->AddBin(5, x11, y11);
    poly->AddBin(6, x12, y12);
    poly->AddBin(8, x13, y13);

    poly->SetLineStyle(1);
    poly->SetLineWidth(10);

    return poly;
 }


///////////////////////////////////////////////////////////////////////////////////////////////////

 //Function generates a vector of TLines to draw the outline for the current binning
 vector<TLine> makeMarsOutline(){
    vector<TLine> tvec(52);

vector<string> cvec = {
  "",
  "Q2 < 2.658 && xB < 0.18",
  "Q2 < 2.658 && xB > 0.18 && xB < 0.21",
  "Q2 < 2.658 && xB > 0.21 && xB < 0.24",
  "Q2 < 2.658 && xB > 0.24 && xB < 0.28",
  "Q2 < 2.658 && xB > 0.28",

  "Q2 > 2.658 && Q2 < 3.625 && xB < 0.24",
  "Q2 > 2.658 && Q2 < 3.625 && xB > 0.24 && xB < 0.28",
  "Q2 > 2.658 && Q2 < 3.625 && xB > 0.28 && xB < 0.34",
  "Q2 > 2.658 && Q2 < 3.625 && xB > 0.34",

  "Q2 > 3.625 && Q2 < 5.12 && xB < 0.34",
  "Q2 > 3.625 && Q2 < 5.12 && xB > 0.34 && xB < 0.40",
  "Q2 > 3.625 && Q2 < 5.12 && xB > 0.40",

  "Q2 > 5.12"
};

    //bin 1
    tvec[0] = TLine(0.126, 2, 0.18, 2.658); //.18
    tvec[1] = TLine(0.18, 2.658, 0.18, 2); //.18
    tvec[2] = TLine(0.18, 2., 0.126, 2);     //.18

    //bin 2
    tvec[3] = TLine(0.18, 2, 0.18, 2.658); //.18
    tvec[4] = TLine(0.18, 2.658, 0.21, 2.658); //.18
    tvec[5] = TLine(0.21, 2.658, 0.21, 2.); 
    tvec[6] = TLine(0.21, 2., 0.18, 2.); //.18

    //bin 3
    tvec[7] = TLine(0.21, 2, 0.21, 2.658); //.18
    tvec[8] = TLine(0.21, 2.658, 0.24, 2.658); //.18
    tvec[9] = TLine(0.24, 2.658, 0.24, 2.); 
    tvec[10] = TLine(0.24, 2., 0.21, 2.); //.18

    //bin 4
    tvec[11] = TLine(0.24, 2., 0.24, 2.658);//.18 
    tvec[12] = TLine(0.24, 2.658, 0.28, 2.658); 
    tvec[13] = TLine(0.28, 2.658, 0.28, 2.); //.18
    tvec[14] = TLine(0.28, 2., 0.24, 2.); //.18

    //bin 5
    tvec[15] = TLine(0.28, 2.658, 0.45, 2.658); 
    tvec[16] = TLine(0.45, 2.658, 0.387826, 2.); 
    tvec[17] = TLine(0.387826, 2., 0.28, 2); 

    //bin 6
    tvec[18] = TLine(0.18, 2.658, 0.24, 3.625);//.18 
    tvec[19] = TLine(0.24, 3.625, 0.24, 2.658); 
    tvec[20] = TLine(0.24, 2.658, 0.18, 2.658); //.18

    //bin 7 
    tvec[21] = TLine(0.24, 2.658, 0.24, 3.625); 
    tvec[22] = TLine(0.24, 3.625, 0.28, 3.625); 
    tvec[23] = TLine(0.28, 3.625, 0.28, 2.658); 
    tvec[24] = TLine(0.28, 2.658, 0.24, 2.658);  

    //bin 8
    tvec[25] = TLine(0.28, 2.658, 0.28, 3.625); 
    tvec[26] = TLine(0.28, 3.625, 0.34, 3.625); 
    tvec[27] = TLine(0.34, 3.625, 0.34, 2.658); 
    tvec[28] = TLine(0.34, 2.658, 0.28, 2.658);  

    //bin 9 
    tvec[29] = TLine(0.34, 2.658, 0.34, 3.625); 
    tvec[30] = TLine(0.34, 3.625, 0.54025, 3.625); 
    tvec[31] = TLine(0.54025, 3.625, 0.45, 2.658);    
    tvec[32] = TLine(0.45, 2.658, 0.34, 2.658);  

    //bin 10
    tvec[33] = TLine(0.24, 3.625, 0.34, 5.12); 
    tvec[34] = TLine(0.34, 5.12, 0.34, 3.625); 
    tvec[35] = TLine(0.34, 3.625, 0.24, 3.625); 

    //bin 11 
    tvec[36] = TLine(0.34, 3.625, 0.34, 5.12); 
    tvec[37] = TLine(0.34, 5.12, 0.40, 5.12); 
    tvec[38] = TLine(0.40, 5.12, 0.40, 3.625);
    tvec[39] = TLine(0.40, 3.625, 0.34, 3.625);

    //bin 12 
    tvec[40] = TLine(0.40, 3.625, 0.40, 5.12); 
    tvec[41] = TLine(0.40, 5.12, 0.6234, 5.12); 
    tvec[42] = TLine(0.6234, 5.12, 0.57, 4.05);    
    tvec[43] = TLine(0.57, 4.05, 0.54025, 3.625);
    tvec[44] = TLine(0.54025, 3.625, 0.40, 3.625);

    //bin 13
    tvec[45] = TLine(0.34, 5.12, 0.677, 10.185); 
    tvec[46] = TLine(0.677, 10.185, 0.7896, 11.351);
    tvec[47] = TLine(0.7896, 11.351, 0.75, 9.52); 
    tvec[48] = TLine(0.75, 9.52, 0.708, 7.42); 
    tvec[49] = TLine(0.708, 7.42, 0.64, 5.4); 
    tvec[50] = TLine(0.64, 5.4, 0.6234, 5.12);
    tvec[51] = TLine(0.6234, 5.12, 0.34, 5.12);

    for (int i = 0; i < tvec.size(); i++){
      tvec[i].SetLineWidth(2);
    }

    return tvec;
 }

