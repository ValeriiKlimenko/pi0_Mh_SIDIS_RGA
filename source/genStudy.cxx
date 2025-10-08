#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <fstream>

#include <boost/range/combine.hpp>

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TError.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TRatioPlot.h"
#include "THStack.h"
#include "THn.h"

#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
//#include "cuts.cxx"



// 01/16/2025 2:35 pm: 
// sampleType::unf IS RESULT OF make_rec_unfold_matricies.cxx some other comments about it can be incorrect
// it is required for REC events 

// Number of THE LAST of x-q2 bins that is used in Dataframe.cxx to assign bins Numbers
// Most of the time it should be used as N_xq2bins + 1 cause it does not take '0' into account
// It is also used in run_unf_matrix for iteration over different bins

#include "binning_params.cxx"
#include "functions_fitting.cxx"
#include "Dataframe.cxx"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;
using namespace std;
using namespace RooFit;

double pi_n = TMath::Pi();

const saveFolder = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/out_filter/gen_test/" 


int countLinesInFile(const std::string& filename) {
    std::ifstream file(filename);
    int lineCount = 0;
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return -1;
    }

    while (std::getline(file, line)) {
        ++lineCount;
    }
    file.close();
    return lineCount;
}

/////////////////////////////////////////////////////////////////

// Change Z-Pt binning




void run_unfold_matrix_new(const int i_batch_files, const int n_files_per_batch, const int iProcedureStep){

///// Setup files range /////////////////////// 
  
  //obtain total number of files
  int max_files = countLinesInFile("path_to_files/nSidis_45_50nA.dat");
  const int n_max_batches = max_files / n_files_per_batch + ( max_files % n_files_per_batch == 0 ? 0 : 1);
  
  const int startFile = n_files_per_batch * i_batch_files;
  if (startFile > max_files)  throw out_of_range("i_batch_files is out of range. There are no files for this batch.")
    
  int endFile = startFile + n_files_per_batch;
  if (endFile > max_files) endFile = max_files;
  if (iProcedureStep < 0 || iProcedureStep > 2) throw out_of_range("iProcedureStep is out of range. 
    The code has two steps: iProcedureStep = 0 for binning generating, 
    iProcedureStep = 1 for filling, and iProcedureStep = 2 is not used");

  TH1::SetDefaultSumw2();
  gROOT->SetBatch(kTRUE);
  
  vector<string> mcRecFiles;
  mcRecFiles = getRecPaths();

  int gstart = 0;
    
  cout <<  " getting RDF: 0%, ";
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame rdf0("h22", mcRecFiles);
  ROOT::RDataFrame rdf_binning("h22", mcRUnFile);
  cout << " RDF is ready (100%)" << endl;
  
// Kinem Functions ////////////////////////////////////////////////////////////////////////////////

  // z bins decreaed by Valerii from 8 to 6 (DO NOT USE IT, IT DOS NOT MATCH OTHER BINNING)
  // it was probably used for some studies.
  auto refz   = new TH1D("refz",   "", 6, 0.2, 0.8);
  auto refp   = new TH1D("refp",   "", 8, 0, 1);
  //There are 9 bins in phi now.
  auto reft12 = new TH1D("reft12", "", 12, 0, 360);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// VALERII BINNING: ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (int i = 0; i < N_xq2bins + 1; i++){
    zpPolyVec[i] = makeTH2PolyZpt2_Valerii(i, zbins_Valerii(i)); 
    zpPolyVec[i]->SetName(("poly_"+ to_string(i)).c_str());
  }

// Cuts that have to capture parameters /////////////////////////////////////////
///// It is better to make parameters global and change lambda -> normal functions

 auto cut_SF=[&]( const double &sf, const double &p, const int &sec, const int strictness){

    // 3.5 nominal (strict = 2)
    double N_sigmas = (5.5 - strictness);
    
    double sf_sigma[3][7] = {{0.01094, 0.01053, 0.01002, 0.01274, 0.01799, 0.01182, 0.01604},
     {0.03428, 0.02436, 0.02993, 0.01008, -0.03193, 0.01386, -0.01597},
     {-0.04474, -0.01564, -0.03043, 0.01254, 0.08689, 0.00483, 0.04774}};

    double sf_mean[3][7] = {{0.20005, 0.20945, 0.1962, 0.23018, 0.22257, 0.1981, 0.22105},
    {0.29526, 0.27354, 0.39302, 0.15601, 0.17654, 0.35734, 0.22785},
    {-0.54065, -0.50649, -0.72252, -0.33457, -0.34202, -0.64477, -0.44011}};

    int sec_for_calc = sec - 1;
    // last is for simulation
    if (isMC) sec_for_calc = 6;

    if (sec_for_calc<0) throw out_of_range("wrong sector in SF cut");

    const double mean = sf_mean[0][sec_for_calc] + sf_mean[1][sec_for_calc] / p  + sf_mean[2][sec_for_calc] / (p*p);
    const double sigma = sf_sigma[0][sec_for_calc] + sf_sigma[1][sec_for_calc] / p  + sf_sigma[2][sec_for_calc] / (p*p);
    
    return sf > mean - sigma * N_sigmas;
};

////////////////////////////////////////
//Dataframe Defin.
////////////////////////////////////////
    
  auto r1 = ROOT::RDF::RNode(rdf0);
  r1 = AddDefine_Kinematics(r1);
  r1 = AddDefine_Kinematics_RecData(r1);
  r1 = AddDefine_Kinematics_RecOnly(r1);

  auto t1 = ROOT::RDF::RNode(rdf_binning);
  t1 = AddDefine_Kinematics_UnfOnly(t1);
  

  //Cut Lists////////////////////////////////////////////////////////////////////////
  // Add PCAL and DC cuts to both main and DIS
  string newCuts       = GetMainCuts(true);//          += "&& dcSecCutEle==1";//Dcsector cut for electron

  //Dis electron cuts
  string rga_fid_cuts  = "y < 0.8 && W > 2 && Q2 > 1 && (e_mom > 2 && e_mom < 8)";
  string dis_mul_cuts = rga_fid_cuts;
  dis_mul_cuts       += "&& y < 0.75 && Q2 > 1";//change q2 > 2 later
  
  //Match cuts
  string match_cut     = "eVal < 0.01 && g1mVal < 0.01 && g2mVal < 0.01";
  string xor_match_cut = "eVal < 0.01 && (g1mVal < 0.01 ^ g2mVal < 0.01)";

  //MC 4 component cuts
  string sig_cut       = "g1mPID==22 && g2mPID==22 && g1mPPID==111 && g2mPPID==111 && g1mPIndex==g2mPIndex";
  string back_cut      = "!(" + sig_cut + ")";
  string tt_cut        = sig_cut;
  string ttf_cut       = "g1mPID==22 && g2mPID==22 && g1mPPID==111 && g2mPPID==111 && g1mPIndex != g2mPIndex";
  string tf_cut        = "(g1mPID==22 ^ g2mPID==22) && (g1mPPID==111 ^ g2mPPID==111)";
  string ff_cut        = "g1mPID!=22 && g2mPID!=22";
  //End cuts/////////////////////////////////////////////////////////////////////////

  //Ancillary Cuts
  //Note that for the unfolding matricies, I am dropping any Q2, x, z, pt2 cuts, but only analyzing the events 
  //that make the eventual cut. However, the explicit xQ2zpt2 bins enforce the cuts with space for underflow and
  //overflow for z bins and overflow bins for pt2. There is also a bin 0 for xQ2 to handle events outside the binning
  //which would be events with x < 0.1 && x > 0.8 and/or events with Q2 < 2.

  string bin_cuts_dis = "Q2>1 && Q2<12";//Accounting for migration
  string bin_cuts_mul = bin_cuts_dis;

////////////////////////////////////////
//CUTS APPLICATION:
////////////////////////////////////////
  
  auto r3 = r1;
  auto r2 = r1;
  r2 = r1.Filter(newCuts.c_str()).Filter(bin_cuts_mul.c_str()).Filter("g1match*g2match>0");
  r3 = r2.Filter(tt_cut.c_str()).Filter("pi0_m > 0.1 && pi0_m < 0.168");


  auto t2 = t1;
  t2 = t1.Filter(bin_cuts_mul.c_str());

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Defining the number of z, pt2, phi trento bins
  //10 z bins : underflow (z < 0.2), 8 bins (0.2 < z < 0.8) , overflow(z > 0.8)
  //9 pt2 bins : 8 bins (0 < pt2 < 1), 1 overflow (pt2 > 1)
  //9 phi trento bins : 9 bins (0 < phi_trento < 360)
  //Note the z and pt2 bins are unique to each xQ2 bin as outlined in zbins and makeTH2PolyZpt2
  // Marshalls:
  //int zpt_bins = 9*10*9+9;
  //Valerii:
  // 10 z bins - > 8 z bins
  // 8+1 pt2 bins -> 13+1 pt2 bins
  // Valerii's hardCoded:
  //int zpt_bins = 9*8*14+9;
  // N_Zbins already includes under and overflow
  int zpt_bins = N_phiTrbins * (N_Zbins) * (N_pTbins + 1) + N_phiTrbins;

  
   bool compute = (iProcedureStep == 0 || iProcedureStep == 0);
   //Bin Center functions
   //The first section makes 2D histograms of all of the phi trento bins by 100 bins in z or pt2.
   //CAN BE DONE IN BATCHES
   if (compute){

    std::vector<std::shared_ptr<TH2D>> hist_z;
    std::vector<std::shared_ptr<TH2D>> hist_pt;

    for (int i = 1; i <= N_xq2bins; ++i) {
        std::string cut_expr = "bin_xBQ2_Valerii == " + std::to_string(i);
    
        hist_z.push_back(
            r3.Filter(cut_expr)
               .Histo2D({Form("r9z_%d", i), "", zpt_bins, 0.5, (double)zpt_bins + 0.5, 100, 0, 1.0},
                        "zpt2phit_8x8x9", "z")
        );
    
        hist_pt.push_back(
            r3.Filter(cut_expr)
               .Histo2D({Form("r9p_%d", i), "", zpt_bins, 0.5, (double)zpt_bins + 0.5, 100, 0, 1.0},
                        "zpt2phit_8x8x9", "pi0_sidis_PT2")
        );
    }
    TString rname = Form("%spi0_TT_bincenter_3d_%d.root", saveFolder.c_str(), startFile);
    auto tfn = TFile::Open(rname.Data(), "recreate");
    for (auto& h : hist_z) h->Write();
    for (auto& h : hist_pt) h->Write();
    tfn->Close();
     
  }

  //Bin Center calculation
  //This function takes the above rootfiles and combines the histograms. From there it finds the bin center for each z-pt2 and z bin, and saves it to a 
  //new root file.
  // SHOULD BE DONE IN ONE GO (NO BATCHES)
  compute = (iProcedureStep == 1 || iProcedureStep == 1);
  //cout<<compute<<" type to int:"<<static_cast<int>(type) << " bf to int:"<<static_cast<int>(sampleType::binning_filling) <<endl;
  if (compute){
    //cout<<"This function takes the above rootfiles and combines the histograms"<<endl;
    //Gathering the files
    //vector<int> num = {0, 30, 60, 90, 120, 150, 180};//, 210};
    
    
    vector<TFile*>fvec(n_max_batches);
    
    for (int i = 0; i < fvec.size(); i++){
      TString rname = Form("%spi0_TT_bincenter_3d_%d.root", saveFolder.c_str(), i_batch_files*i);
      fvec[i] = new TFile(rname.Data());
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Redone to make it shorter and more efficient: ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //Histogram vector initialization and adding
    //Rebining by 9 integrates over phi trento and if done another time it integrates over pt2
    vector<TH2D*> hvz(N_xq2bins + 1), hvp(N_xq2bins + 1), hvz2(N_xq2bins + 1);
    vector<TH1D*> hvz_out(N_xq2bins + 1), hvp_out(N_xq2bins + 1), hvz2_out(N_xq2bins + 1);

    for (int i = 1; i <= N_xq2bins; i++) {
      string idx = to_string(i);
      hvz[i] = (TH2D*)fvec[0]->Get(("r9z_" + idx).c_str())->Clone();
      hvp[i] = (TH2D*)fvec[0]->Get(("r9p_" + idx).c_str())->Clone();

      // the N bins is the same because it is zpt2phit_8x8x9 for both histograms
      int n_hv_integral = hvz[i]->GetNbinsX() / N_phiTrbins;
    
      // TO DO: It would be better to use dimensions from rebinned histograms. It wound eliminate all posiible +-1 rounding issues.
      hvz_out[i]  = new TH1D(("r9z_cen_" + idx).c_str(), "", n_hv_integral, 0.5, n_hv_integral + 0.5);
      hvp_out[i]  = new TH1D(("r9p_cen_" + idx).c_str(), "", n_hv_integral, 0.5, n_hv_integral + 0.5);
      hvz2_out[i] = new TH1D(("r9z_int_cen_" + idx).c_str(), "", n_hv_integral - 1, 0.5, n_hv_integral - 1 + 0.5);
    }

    //Adding the histograms from all the files for every xQ2 bin
    for (size_t i = 1; i < fvec.size(); ++i) {
      for (int j = 1; j <= N_xq2bins; ++j) {
        hvz[j]->Add((TH2D*)fvec[i]->Get(("r9z_" + to_string(j)).c_str()));
        hvp[j]->Add((TH2D*)fvec[i]->Get(("r9p_" + to_string(j)).c_str()));
      }
    }

    // Rebinning the histograms by 9, which integrates over phi_trento
    for (int i = 1; i <= N_xq2bins; ++i) {
      hvz[i]->RebinX(N_phiTrbins);
      hvp[i]->RebinX(N_phiTrbins);

      // the same for hvz and hvp
      int n_rebinned_bins = hvz[i]->GetNbinsX();
    
      for (int j = 0; j < n_rebinned_bins; ++j) {
        TH1D lhz("", "", 100, 0, 1);
        TH1D lhp("", "", 100, 0, 1);
    
        for (int k = 1; k <= 100; ++k) {
          lhz.SetBinContent(k, hvz[i]->GetBinContent(j + 1, k));
          lhz.SetBinError(k,   hvz[i]->GetBinError(j + 1, k));
          lhp.SetBinContent(k, hvp[i]->GetBinContent(j + 1, k));
          lhp.SetBinError(k,   hvp[i]->GetBinError(j + 1, k));
        }
    
        hvz_out[i]->SetBinContent(j + 1, lhz.GetMean());
        hvp_out[i]->SetBinContent(j + 1, lhp.GetMean());
    
        if (j < n_rebinned_bins - 1 && j <= N_pTbins + 1) {
          TH1D lhz2("", "", 100, 0, 1);
          for (int k = 1; k <= 100; ++k) {
            lhz2.SetBinContent(k, hvz[i]->GetBinContent(j + 2, k));
            lhz2.SetBinError(k,   hvz[i]->GetBinError(j + 2, k));
          }
          hvz2_out[i]->SetBinContent(j + 1, lhz2.GetMean());
        }
      }
    }

    TString rname = Form("%spi0_TT_bincenter_results.root", saveFolder.c_str());
    auto tfout = new TFile(rname.Data(), "recreate");
    for (int i = 1; i <= N_xq2bins; ++i) {
      hvz_out[i]->Write();
      hvz2_out[i]->Write();
      hvp_out[i]->Write();
    }
    
 }

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //Making 1st step matricies
 //This generates histograms of the binning for the gen, rec, data, and dis
 //For data, the histogram is 2D, with the flattened 3D binning on the x axis and the mgg mass distribution on the y axis
 //For rec, the histogram is 3D, with the flattended 3D binning on the x axis, the flattend 3D binning of the matched events
 //on the y axis, and the mgg distribution on teh z axos
 //For gen, the histogram is simply a 1D histigram of the flattend 3D binning
 //For dis, the histogrm is a 1D histogram of the xQ2 bins


  //Rec
  // DOES NOT NEED RESULTS OF iProcedureStep = 0,1
  compute = (iProcedureStep == 2 || iProcedureStep == 2);
  if (compute == 1){

    TString rname = Form("%spi0_rec_v3_unfold_%d.root", filePath.c_str(), startFile);
    auto tfn = new TFile(rname.Data(), "recreate");
    for (int i = 1; i < N_xq2bins + 1; i++){
      string cut = "bin_xBQ2_Valerii ==" + to_string(i);
      auto r9  = t2.Filter(cut.c_str()).Histo3D({("r9_" + to_string(i)).c_str(),  "", zpt_bins, 1, (double)zpt_bins+1, zpt_bins, 1, (double)zpt_bins+1, N_pi0mm_bins, min_pi0mass, max_pi0mass}, "zpt2phit_8x8x9", "zpt2phit_8x8x9m", "pi0_m");
      r9->Write();
    }
    tfn->Close();
    cout << rname << endl;
  }




  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Fiting pi0 distributions
  //These routines fit the pi0 distributions within the data and rec histograms, and then generate new histograms with the fit results
  //So for data the 2D histogram in bins vs mgg becomes a 1D histogram in just bins with the bin content being the number of extracted pi0s
  //For rec, the 3D histogram of bins vs matched matched bins vs mgg becomes bins vs matched bins with the bin content being the number of extracted pi0s
  //The errors from  the fits are also propagated and become the errors on the bin content
  //For reference the *vec[] are the full z-pt2-phitrento bins and the *vec2[] are the phi_trento integrated bins, so z-pt2 bins
  //Histogram nomenclature:
  //Data
  //d9_*      are histograms for the full 3D binning, fit
  //dpi0_*    are histograms for the 2D binning, phitrento integrated bins, fit
  //dsum_*    are unfit histograms for the full 3D binning, that have been projected down to 1D histograms, unfit
  //d9_2d_*   are histograms are 2D histograms with bins vs mgg, unfit
  //Rec
  //r9_*      are histograms for the full 3D binning, fit
  //rpi0_*    are histograms for the 2D binning, phitrento integrated bins, fit
  //rsum_*    are unfit histograms for the full 3D binning, that have been projected down to 1D histograms, unfit
  //r9_2d_*   are histograms are 2D histograms with bins vs mgg, unfit
  //rpi0_1D_* are histograms for the full 3D binning, matching integrated, then fit for a resultant 1D histogram
  //Gen
  //g9_*      are histograms for the full 3D binning, unfit
  //gpi0_*    are histograms for the 2D binning, phitrento integrated bins, unfit

  //Data
  compute =  (type == sampleType::data && second_fitted);
  if (compute == 1){

    vector<TH2D*> dvec(N_xq2bins);
    vector<TH2D*> dvec2(N_xq2bins);
    for (int i = 1; i < N_xq2bins + 1; i++){
      TString dname = Form("%spi0_data_v3_unfold_0.root", filePath.c_str());
      //cout<<dname<<endl;
      auto tf0 = new TFile(dname.Data());
      dvec[i-1] = (TH2D*)((TH2D*)tf0->Get(("d9_" + to_string(i)).c_str()))->Clone();
      dvec2[i-1] = (TH2D*)dvec[i-1]->RebinX(N_phiTrbins, ("dpi0_" + to_string(i)).c_str());
    }

    for (int i = 50; i < 156; i+=50){
      TString dname = Form("%spi0_data_v3_unfold_%d.root", filePath.c_str(), i);
      auto tf = new TFile(dname.Data());
      for (int j = 1; j < N_xq2bins + 1; j++){
        dvec[j-1]->Add((TH2D*)tf->Get(("d9_" + to_string(j)).c_str()));
        auto h = (TH2D*)((TH2D*)tf->Get(("d9_" + to_string(j)).c_str()))->Clone();
        // integration over phi:
        h->RebinX(N_phiTrbins);
        dvec2[j-1]->Add(h);
      }
    }

    TString name = Form("%spi0_data_v3_unfold_phit9_1D.root", filePath.c_str());
    auto tfn = new TFile(name.Data(), "recreate");
    
    for (int i = 1; i < N_xq2bins + 1; i++){
    //for (int i = 1; i < 2; i++){
      auto h = extractPi0(dvec[i-1]);
      h->SetName(("d9_" + to_string(i)).c_str());
      h->Write();
      // dvec2 is an integrated dvec
      // Input is a vector of histograms, where every vector's element represents xQ2 bin, the element itself is 2D histogram where
      // X axis is linerlized binning for three variables (int zpt_bins = 9*8*14+9;) z, pt and phi,  pi0_m is Y variable
      // "h->RebinX(9)" integrates over phi so we left with X axis that represents Z and Pt
      auto h2 = extractPi0(dvec2[i-1], n_sigma_pi0mass, true);
      h2->SetName(("dpi0_" + to_string(i)).c_str());
      h2->Write();
      auto h3 = dvec2[i-1]->ProjectionX();
      h3->SetName(("dsum_" + to_string(i)).c_str());
      h3->Write();
      auto h4 = (TH2D*)dvec[i-1]->Clone();
      h4->SetName(("d9_2d_" + to_string(i)).c_str());
      h4->Write();
      
    }
    
  }

  //Rec
  compute = (type == sampleType::rec && second_fitted);
  if (compute == 1){

    vector<TH3D*> rvec(N_xq2bins);
    vector<TH3D*> rvec2(N_xq2bins);
    for (int i = 1; i < N_xq2bins + 1; i++){
      TString rname = Form("%spi0_rec_v3_unfold_0.root", filePath.c_str());
      auto tf0 = new TFile(rname.Data());
      rvec[i-1] = (TH3D*)((TH3D*)tf0->Get(("r9_" + to_string(i)).c_str()))->Clone();
      rvec2[i-1] = (TH3D*)rvec[i-1]->RebinX(N_phiTrbins, ("rpi0_" + to_string(i)).c_str());
      // integration over generated trento phi
      rvec2[i-1]->RebinY(N_phiTrbins);

    }

    for (int i = 30; i < maxMCfiles; i+=30){
      TString rname = Form("%spi0_rec_v3_unfold_%d.root", filePath.c_str(), i);
      auto tf = new TFile(rname.Data());
      for (int j = 1; j < N_xq2bins + 1; j++){
        rvec[j-1]->Add((TH3D*)tf->Get(("r9_" + to_string(j)).c_str()));
        auto h = (TH3D*)((TH3D*)tf->Get(("r9_" + to_string(j)).c_str()))->Clone();
        h->RebinX(N_phiTrbins);
        // integration over generated trento phi
        h->RebinY(N_phiTrbins);
        rvec2[j-1]->Add(h);
      }
    }

    TString name = Form("%spi0_rec_v3_unfold_phit9_1D.root", filePath.c_str());
    auto tfn = new TFile(name.Data(), "recreate");
    //auto tfn = new TFile(name.Data(), "update");
    
    for (int i = 1; i < N_xq2bins + 1; i++){
      auto h = extractPi02D(rvec[i-1], n_sigma_pi0mass, true);
      cout << "RVEC EXTRACTED"<<endl;
      h->SetName(("r9_" + to_string(i)).c_str());
      h->Write();
      auto h2 = extractPi02D(rvec2[i-1], n_sigma_pi0mass, true);
      cout << "RVEC2 EXTRACTED"<<endl;
      h2->SetName(("rpi0_" + to_string(i)).c_str());
      h2->Write();
      auto h3 = rvec2[i-1]->Project3D("xy");
      h3->SetName(("rsum_" + to_string(i)).c_str());
      h3->Write();
      auto h4 = (TH3D*)rvec[i-1]->Clone(("r9_3d_" + to_string(i)).c_str());
      auto h5 = (TH2D*)rvec[i-1]->Project3D("zx");
      h5->SetName(("r9_2d_" + to_string(i)).c_str());
      h4->Write();
      h5->Write();
      auto h6 = (TH2D*)h5->Clone();
      h6->RebinX(N_phiTrbins);
      auto h7 = extractPi0(h6, n_sigma_pi0mass, true);
      cout << "RVEC - ZX proj EXTRACTED"<<endl;
      h7->SetName(("rpi0_1D_" + to_string(i)).c_str());
      h7->Write();

      printf("Done # %d bin out of %d bins\n",i,N_xq2bins + 1);
    }
  }

  //Gen
  compute = (type == sampleType::gen && second_fitted);
  if (compute == 1){
    vector<TH1D*> gvec(N_xq2bins);
    vector<TH1D*> gvec2(N_xq2bins);
    for (int i = 1; i < N_xq2bins + 1; i++){
      TString gname = Form("%spi0_gen_v3_unfold_0.root", filePath.c_str());
      auto tf0 = new TFile(gname.Data());
      gvec[i-1] = (TH1D*)((TH1D*)tf0->Get(("g9_" + to_string(i)).c_str()))->Clone();
      gvec2[i-1] = (TH1D*)gvec[i-1]->RebinX(N_phiTrbins, ("gpi0_" + to_string(i)).c_str());
    }

    for (int i = 30; i < maxMCfiles; i+=30){
      TString gname = Form("%spi0_gen_v3_unfold_%d.root", filePath.c_str(), i);
      auto tf = new TFile(gname.Data());
      for (int j = 1; j < N_xq2bins + 1; j++){
        gvec[j-1]->Add((TH1D*)tf->Get(("g9_" + to_string(j)).c_str()));
        auto h = (TH1D*)((TH1D*)tf->Get(("g9_" + to_string(j)).c_str()))->Clone();
        h->RebinX(N_phiTrbins);
        gvec2[j-1]->Add(h);
      }
    }

    TString name = Form("%spi0_gen_v3_unfold_phit9_1D.root", filePath.c_str());
    auto tfn = new TFile(name.Data(), "recreate");
    for (int i = 1; i < N_xq2bins + 1; i++){
      gvec[i-1]->Write();
      gvec2[i-1]->SetName(("gpi0_" + to_string(i)).c_str());
      gvec2[i-1]->Write();
    }
  }

}//file end
