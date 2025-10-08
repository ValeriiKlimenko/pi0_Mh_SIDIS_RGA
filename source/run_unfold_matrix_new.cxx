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

#include "momentum_corrections.cxx"

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

/////////////////////////////////////////////////////////////////

// Change Z-Pt binning

void run_unfold_matrix_new(int sPoint, int type_of_compute_INT, bool second_fitted = false){

  if (second_fitted){
    // Error Ignore Level Set
    gErrorIgnoreLevel = kFatal;
    cout << "Error Ignore Level Set to kFatal because of fit output" <<endl;
  }

  if (type_of_compute_INT > 5 || type_of_compute_INT < 0) 
    throw out_of_range("Wrong type of files in sampleType");
  
  sampleType type = static_cast<sampleType>(type_of_compute_INT);
  check_input_run_unf(type, sPoint);

  // Parameters: /////
  // This one is for SF cut:
  bool isMC = (type == sampleType::rec || type == sampleType::gen || type == sampleType::binning_initial) ? true : false;
  
  //Check the sPoint and data type
  //Data files are being initialized even when the MC calculations are running and the data is not being used. So, this is a fix for now:
  const int sPointData = (sPoint > 150) ? 150 : sPoint;
  cout<<"Start File: "<<sPoint<<endl;

  //files limit
  const int maxMCfiles = 210;
  const int maxDataFiles = 171;
  TH1::SetDefaultSumw2();

  //Data : 174 files, Rec and Gen 220 files for now
  vector <string> dataFile, mcRecFile, mcGenFile, mcRUnFile, disFile, tdata, tmcR, tmcG;

/////////////////////////// Data files reading /////////////////////
  // only one vector out od four will contatin information (type will determine which one)
  dataFile = GetDATA_filePaths(type);
  disFile = GetDATA_filePaths(type);

  mcRecFile = GetMC_filePaths(type);
  mcGenFile = GetMC_filePaths(type);
/////////////////////////// Subset of Files: //////////////////////////  

  //For running the files, the files are broken into groups of files
  //groups of 50 for data and dis, and groups of 30 for gen and rec. 
  //Later functions open all of the files to sum them and/or fit distributions;

 int dstart = 0, sstart = 0, gstart = 0, rstart = 0, mcint = 0;


 cout <<  " get subset " << endl;
  
 GetSubsetOfDataFiles(dataFile, mcRecFile, mcGenFile, mcRUnFile, disFile,
   sPointData, sPoint, maxDataFiles, maxMCfiles,
   // starting point for each file type
   dstart, sstart, gstart, rstart, mcint
  );

    cout <<  " get RDF " << endl;
////////// Data Frames initialization (Only one is actually filled (determined by "type") /////////////////
  ROOT::EnableImplicitMT();

  // in order to avoid the sitation when there is no files to open the ttree, just one RDF is used for all types:
  auto filePathsRDF=[&](const sampleType type){
    if (type == sampleType::rec || type == sampleType::binning_initial ) return mcRecFile;
    if (type == sampleType::gen) return mcGenFile;
    // the files are created after inital running for the programm "make binning", it is required for REC files as a replacement of mcRecFile(?)
    // I do not know for sure so I will open both mcRecFile and mcRUnFile in case of sampleType::REC
    if (type == sampleType::unf) return mcRUnFile;
    if (type == sampleType::data) return dataFile;
    if (type == sampleType::dis) return disFile;
  
    return vector<string>();
  };

  
  ROOT::RDataFrame rdf0("h22", filePathsRDF(type));
  // Hot fix to have "mcRUnFile" for REC type. It will be initilliazed always so it will be avaliable for ::REC.
  ROOT::RDataFrame rdf_binning("h22", mcRUnFile);
  
/*
  ROOT::RDataFrame r0("h22", mcRecFile);
  ROOT::RDataFrame g0("h22", mcGenFile);
  ROOT::RDataFrame t0("h22", mcRUnFile);
  ROOT::RDataFrame d0("h22", dataFile);
  ROOT::RDataFrame s0("h22", disFile);
*/
  
  // Kinem Functions ////////////////////////////////////////////////////////////////////////////////
  cout <<  " after RDF " << endl; 
  
/*
 //Generates xQ2 outlines
 auto out_mars = makeMarsOutline();

  //Returns bin number in the xB vs. Q^2 binning developed by Marshall
  // see Table 4.2 in clas12 Sidis analysis note
  auto bin_xBQ2_Mars=[&](const double &x, const double &Q2){
    return poly_mars->FindBin(x, Q2);
 };
*/

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

/*
  auto bin_xBQ2_Valerii=[&](const double &x, const double &Q2){
    return poly_valerii->FindBin(x, Q2);
 };
*/
  
/////////////////////////////////////////////////////
////// FOR WHATEVER REASON IS DONE IN UNIFORM BINNING:  
////// so it is a different binning ///////////////// 
/////////////////////////////////////////////////////

/*  
//Returns zpt2 bin number
 auto zpt2_8x8=[&](double &z, double &pt2){
  // 8 equal bins in pt2, //6 equal bins in z 0.2 - 0.8 (Decreased from 8 by Valerii)
  // 10 -> 8 cause decrease in number of Z bins (from 8 to 6 + 2 under/overflow)
  return 8*refp->FindBin(pt2) + refz->FindBin(z);
};
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END VALERII BINNING: ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
//OLD Returns zpt2 bin number
 auto zpt2phit_8x8x9=[&](double &z, double &pt2, double &phit){
     return 10*11*refp->FindBin(pt2) + 11*refz->FindBin(z) + reft9->FindBin(phit);
};
*/

/// SAME PROBLEM, ALL ARE UNIFORM WHILE BOTH MARSHALL AND I DID NOT USE UNIFORM:
  
/*
//Returns xq2zpt2 bin number
 auto xq2zpt2_13x8x8=[&](double &x, double &Q2, double &z, double &pt2){
     int mars_bin = poly_mars->FindBin(x, Q2);
     if (mars_bin < 1){mars_bin = 0;}
     return 10*10*mars_bin + 10*refp->FindBin(pt2) + refz->FindBin(z);

};

//Returns xq2zpt2phit bin number 12 bins
 auto xq2zpt2phit_13x8x8x12=[&](double &x, double &Q2, double &z, double &pt2, double &phit){
     int mars_bin = poly_mars->FindBin(x, Q2);
     if (mars_bin < 1){mars_bin = 0;}

     return 14*10*10*mars_bin + 10*14*refp->FindBin(pt2) + 14*refz->FindBin(z) + reft12->FindBin(phit);
};

//Returns xq2zpt2phit bin number 9 bins
 auto xq2zpt2phit_13x8x8x9=[&](double &x, double &Q2, double &z, double &pt2, double &phit){
     int mars_bin = poly_mars->FindBin(x, Q2);
     if (mars_bin < 1){mars_bin = 0;}

     return 11*10*10*mars_bin + 10*11*refp->FindBin(pt2) + 11*refz->FindBin(z) + reft9->FindBin(phit);

};
*/

// Cuts that have to capture parameters /////////////////////////////////////////
///// It is better to make parameters global in change lmbda -> normal functions

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


//Dataframe Defin.
////////////////////////////////////////

  //data
  auto d1 = ROOT::RDF::RNode(rdf0);

  if (type == sampleType::data){
    isMC = false;
    d1 = AddDefine_Kinematics(d1);
    d1 = AddDefine_Kinematics_RecData(d1);
  }
 
  //Recon MonteCarlo DataFrame
  //mc rec
  auto r1 = ROOT::RDF::RNode(rdf0);
  if (type == sampleType::rec || type == sampleType::binning_initial){
    
    isMC = true;
    cout<< "File exmp: "<<mcRecFile[0]<<endl;
    r1 = AddDefine_Kinematics(r1);
    r1 = AddDefine_Kinematics_RecData(r1);
    r1 = AddDefine_Kinematics_RecOnly(r1);
  }
 
  //Recon MonteCarlo DataFrame
  //recrec RESULT OF INITIAL PASS OF THIS PROGRAM (making binning)
  auto t1 = ROOT::RDF::RNode(rdf_binning);
  if (type == sampleType::unf || type == sampleType::rec){
    isMC = true;
    // it should have binning information only  so remove kin.
    //t1 = AddDefine_Kinematics(t1);
    t1 = AddDefine_Kinematics_UnfOnly(t1);
  }
  //Gen Data frame
  //gengen

  auto g1 = ROOT::RDF::RNode(rdf0);
  if (type == sampleType::gen){
    isMC = true;
    g1 = AddDefine_Kinematics(g1);
    g1 = AddDefine_Kinematics_GenOnly(g1);  
  }

  //Dis Data frame
  //disdis

  auto s1 = ROOT::RDF::RNode(rdf0);
  if (type == sampleType::dis){
    isMC = false;
    s1 = AddDefine_Kinematics(s1);
  }
  

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
  //bin_cuts_dis       += "&& xB > 0.1 && xB < 0.8";//Accounting for migration
  string bin_cuts_mul = bin_cuts_dis;

  //bin_cuts_mul       += "&& z > 0.2 && z < 0.8";//0.1-0.8
  //bin_cuts_mul       += "&& pi0_sidis_PT2 > 0 && pi0_sidis_PT2 < 1";

  auto d3 = d1;
  auto d2 = d1;

  //CUTS APPLICATION:
  
  if (type == sampleType::data){
    d2 = d1.Filter(newCuts.c_str()).Filter(bin_cuts_mul.c_str());
    d3 = d2.Filter("pi0_m > 0.1 && pi0_m < 0.168");
  }  

  auto r3 = r1;
  auto r2 = r1;
  if (type == sampleType::rec || type == sampleType::binning_initial){
    r2 = r1.Filter(newCuts.c_str()).Filter(bin_cuts_mul.c_str()).Filter("g1match*g2match>0");
    r3 = r2.Filter(tt_cut.c_str()).Filter("pi0_m > 0.1 && pi0_m < 0.168");
  }

  // came form make_rec_unfold_matricies.cxx,
  // all the cuts are specified there so here we have Q2 cut only
  // It does hurt performance only to do the cuts there. 
  // However, the cuts need to be sinchronized between this file (run_unfold_matrix_new.cxx) and "make_rec_unfold_matricies.cxx"
  
  auto t2 = t1;
  if (type == sampleType::rec || type == sampleType::unf){
    t2 = t1.Filter(bin_cuts_mul.c_str());
  }
  
  auto s2 = s1;
  // has no cut whatsoever except kinematics
  if (type == sampleType::dis){
    s2 = s1.Filter(dis_mul_cuts.c_str()).Filter(bin_cuts_dis.c_str());
  }

  auto g2 = g1;
  // Valerii: goodRecElec is check on FD, Q2, W and Vz
  if (type == sampleType::gen){
    g2 = g1.Filter("goodRecElec==1").Filter(dis_mul_cuts.c_str()).Filter(bin_cuts_mul.c_str()).Filter("Mx>1.5");
  }

  bool compute = 0;
  gROOT->SetBatch(kTRUE);

  //Testing THnSparseD 
  //auto hj = new TH2D("", "", 100, 0, 1,100, 0, 1);
  //int nbins[3] = {1221, 1221, 20};
  //double mins[3] = {0,0,0};
  //double maxs[3] = {1221, 1221, 20};
  //auto hj2 = new THnSparseD("", "", 3, nbins, mins, maxs);//,{"xq2zpt2phit_13x8x8x9", "xq2zpt2phit_13x8x8x9m", "pi0_m"});
  //r2.Foreach([hj2](int v0, int v1, double v2){hj2->Fill(v0, v1, v2);},{"zpt2phit_8x8x9", "zpt2phit_8x8x9m", "pi0_m"});

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

  
   compute = (type == sampleType::binning_initial);
   //Bin Center functions
   //The first section makes 2D histograms of all of the phi trento bins by 100 bins in z or pt2.
   //This runs for groups of 30 rootfiles.
   //cout<<compute<<" type to int:"<<static_cast<int>(type) << " bf to int:"<<static_cast<int>(sampleType::binning_initial) <<endl;
   if (compute){
    //cout<<"rstart:"<<rstart<<endl;
    TString rname = Form("%spi0_TT_bincenter_3d_%d.root", filePath.c_str(), rstart);
    auto tfn = new TFile(rname.Data(), "recreate");

    // cause of split mow it is 1:16
    for (int i = 1; i < N_xq2bins + 1; i++){
    //for (int i = 1; i < 14; i++){
      //string cut = "bin_xBQ2_Mars==" + to_string(i);
      string cut = "bin_xBQ2_Valerii==" + to_string(i);

      // +1 -> 0.5 to have bin number in the center of the bin
      auto rz  = r3.Filter(cut.c_str()).Histo2D({("r9z_" + to_string(i)).c_str(),  "", zpt_bins, 0.5, (double)zpt_bins+0.5, 100, 0, 1.0}, "zpt2phit_8x8x9", "z");
      rz->Write();
      auto rp  = r3.Filter(cut.c_str()).Histo2D({("r9p_" + to_string(i)).c_str(),  "", zpt_bins, 0.5, (double)zpt_bins+0.5, 100, 0, 1.0}, "zpt2phit_8x8x9", "pi0_sidis_PT2");
      rp->Write();
    }
    tfn->Close();
    cout << rname << endl;
  }

  //Bin Center calculation
  //This function takes the above rootfiles and combines the histograms. From there it finds the bin center for each z-pt2 and z bin, and saves it to a 
  //new root file.
  compute = (type == sampleType::binning_filling);
  //cout<<compute<<" type to int:"<<static_cast<int>(type) << " bf to int:"<<static_cast<int>(sampleType::binning_filling) <<endl;
  if (compute){
    //cout<<"This function takes the above rootfiles and combines the histograms"<<endl;
    //Gathering the files
    vector<int> num = {0, 30, 60, 90, 120, 150, 180};//, 210};
    vector<TFile*>fvec(num.size());
    int ct = 0;
    for (int i = 0; i < num.size(); i++){
      TString rname = Form("%spi0_TT_bincenter_3d_%d.root", filePath.c_str(), num[i]);
      fvec[i] = new TFile(rname.Data());
    }

    //Histogram vector initialization and adding
    //Rebining by 9 integrates over phi trento and if done another time it integrates over pt2
    vector<TH2D*> hvz(N_xq2bins + 1), hvp(N_xq2bins + 1), hvz2(N_xq2bins + 1);
    vector<TH1D*> hvz_out(N_xq2bins + 1), hvp_out(N_xq2bins + 1), hvz2_out(N_xq2bins + 1);

    //Initializing the vectors
    for (int i = 1; i < N_xq2bins + 1; i++){
      //cout<<"i:"<<i<<' '<<fvec[0]<<endl;
      hvz[i] = (TH2D*)((TH2D*)fvec[0]->Get(("r9z_" + to_string(i)).c_str()))->Clone();
      hvp[i] = (TH2D*)((TH2D*)fvec[0]->Get(("r9p_" + to_string(i)).c_str()))->Clone();

      
    
      // Histograms that will be integrated over phi:
      int n_hvz_integral = int( hvz[i]->GetNbinsX() / N_phiTrbins);
      int n_hvp_integral = int( hvp[i]->GetNbinsX() / N_phiTrbins);
      //int n_hvz2_integral = int( hvz2[i]->GetNbinsX() / N_phiTrbins);

      // TO DO: It would be better to use dimensions from rebinned histograms. It wound eliminate all posiible +-1 rounding issues.
      hvz_out[i]  = new TH1D(("r9z_cen_" + to_string(i)).c_str(),"", n_hvz_integral, 0.5, n_hvz_integral + 0.5); 
      hvp_out[i]  = new TH1D(("r9p_cen_" + to_string(i)).c_str(),"",  n_hvp_integral, 0.5, n_hvp_integral + 0.5); 
      //cout<<zpt_bins<<' ' <<hvz[i]->GetNbinsX() << ' '<<n_hvz_integral<<endl;
      hvz2_out[i] = new TH1D(("r9z_int_cen_" + to_string(i)).c_str(),"", n_hvz_integral - 1, 0.5, n_hvz_integral - 1 + 0.5); 

      
    }

    //Adding the histograms from all the files for every xQ2 bin
    for (int i = 1; i < fvec.size(); i++){
      for (int j = 1; j < N_xq2bins + 1; j++){
        hvz[j]->Add((TH2D*)fvec[i]->Get(("r9z_" + to_string(j)).c_str()));
        hvp[j]->Add((TH2D*)fvec[i]->Get(("r9p_" + to_string(j)).c_str()));
      }
    }

    // Rebinning the histograms by 9, which integrates over phi_trento
    for (int i = 1; i < N_xq2bins + 1; i++){
      hvz2[i] = (TH2D*)hvz[i]->Clone(("r9z_int_cen_loc_" + to_string(i)).c_str());

      // TH1F *hnew = dynamic_cast<TH1F*>(h1->Rebin(5,"hnew")); 
      int n_hvz_integral = int( hvz[i]->GetNbinsX() / N_phiTrbins);

      // Integration in place:
      hvz[i]->RebinX(N_phiTrbins);
      hvp[i]->RebinX(N_phiTrbins);
      hvz2[i]->RebinX(N_phiTrbins);
      

      auto h = new TH2D("","", n_hvz_integral - 1, 0.5, n_hvz_integral - 1 + 0.5, 100, 0, 1);
      // Even if number of incorrect we will miss just the last bin
      for (int j = 0; j < n_hvz_integral - 1; j++){
        for (int k = 0; k < 100; k++){
          h->SetBinContent(j+1, k+1, hvz[i]->GetBinContent(j+2, k+1));
          h->SetBinError(j+1, k+1, hvz[i]->GetBinError(j+2, k+1));
        }
      }
      //Integral over pT:
      h->RebinX(N_pTbins + 1);
      hvz2[i] = (TH2D*)h->Clone(("r9z_int_cen_loc_" + to_string(i)).c_str());
    }

    // Finding the means of the histogram bins and writing them to the root file 
    TString rname = Form("%spi0_TT_bincenter_results.root", filePath.c_str());
    auto tfout = new TFile(rname.Data(), "recreate");
    for (int i = 1; i < N_xq2bins + 1; i++){
      for (int j = 0; j < hvz[i]->GetNbinsX(); j++){
        auto lhz = new TH1D("","", 100, 0, 1);
        auto lhp = new TH1D("","", 100, 0, 1);
        auto lhz2 = new TH1D("","", 100, 0, 1);
        for (int k = 0; k < 100; k++){
          lhz->SetBinContent(k+1, hvz[i]->GetBinContent(j+1, k+1));
          lhz->SetBinError(k+1,   hvz[i]->GetBinError(j+1, k+1));
          lhp->SetBinContent(k+1, hvp[i]->GetBinContent(j+1, k+1));
          lhp->SetBinError(k+1,   hvp[i]->GetBinError(j+1, k+1));
        }
        hvz_out[i]->SetBinContent(j+1,lhz->GetMean());
        hvp_out[i]->SetBinContent(j+1,lhp->GetMean());

        // Valerii: hvz2 was integrated ove pt so only Z bins are left for (continue j > 9).
        //if (j > 9){continue;}
        if (j > N_pTbins + 1){continue;}
        for (int k = 0; k < 100; k++){
          lhz2->SetBinContent(k+1, hvz2[i]->GetBinContent(j+1, k+1));
          lhz2->SetBinError(k+1,   hvz2[i]->GetBinError(j+1, k+1));
        }
        hvz2_out[i]->SetBinContent(j+1,lhz2->GetMean());
      }
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
  
  //Data


  compute = (type == sampleType::data && !second_fitted);
  if (compute == 1){

    TString dname = Form("%spi0_data_v3_unfold_%d.root", filePath.c_str(), dstart);
    cout<<"saveFileName:"<<dname.Data()<<endl;
    auto tfn = new TFile(dname.Data(), "recreate");
    for (int i = 1; i < N_xq2bins + 1; i++){
    //for (int i = 1; i < 14; i++){
      //string cut = "bin_xBQ2_Mars==" + to_string(i);
      string cut = "bin_xBQ2_Valerii==" + to_string(i);
      auto d9  = d2.Filter(cut.c_str()).Histo2D({("d9_" + to_string(i)).c_str(),  "", zpt_bins, 1, (double)zpt_bins+1, N_pi0mm_bins, min_pi0mass, max_pi0mass}, "zpt2phit_8x8x9", "pi0_m");
      d9->Write();
    }
    tfn->Close();
    cout << dname << endl;
  }

  //Gen
  compute = (type == sampleType::gen && !second_fitted);
  if (compute == 1){

    TString gname = Form("%spi0_gen_v3_unfold_%d.root", filePath.c_str(), gstart);
    auto tfn = new TFile(gname.Data(), "recreate");
    for (int i = 1; i < N_xq2bins + 1; i++){
    //for (int i = 1; i < 14; i++){
      //string cut = "bin_xBQ2_Mars==" + to_string(i);
      string cut = "bin_xBQ2_Valerii==" + to_string(i);
      auto g9  = g2.Filter(cut.c_str()).Histo1D({("g9_" + to_string(i)).c_str(),  "", zpt_bins, 1, (double)zpt_bins+1}, "zpt2phit_8x8x9");
      g9->Write();
    }
    tfn->Close();
    cout << gname << endl;
  }

  //Rec
  compute = (type == sampleType::rec && !second_fitted);
  if (compute == 1){

    TString rname = Form("%spi0_rec_v3_unfold_%d.root", filePath.c_str(), mcint);
    auto tfn = new TFile(rname.Data(), "recreate");
    for (int i = 1; i < N_xq2bins + 1; i++){
    //for (int i = 1; i < 14; i++){
      //string cut = "bin_xBQ2_Mars==" + to_string(i);
      string cut = "bin_xBQ2_Valerii ==" + to_string(i);
      //string cut = "1 == 1";// + to_string(i);
      auto r9  = t2.Filter(cut.c_str()).Histo3D({("r9_" + to_string(i)).c_str(),  "", zpt_bins, 1, (double)zpt_bins+1, zpt_bins, 1, (double)zpt_bins+1, N_pi0mm_bins, min_pi0mass, max_pi0mass}, "zpt2phit_8x8x9", "zpt2phit_8x8x9m", "pi0_m");
      r9->Write();
    }
    tfn->Close();
    cout << rname << endl;
  }


  //DIS
  compute = (type == sampleType::dis && !second_fitted);
  if (compute == 1){

    TString sname = Form("%sdis_v2_xQ2_stat_%d.root", filePath.c_str(), sstart);
    auto tfn = new TFile(sname.Data(), "recreate");
    //auto s9  = s2.Histo1D({"dis", "", 14, 0, 14}, "bin_xBQ2_Mars");
    //Valerii: check if it is correct scheme
    auto s9  = s2.Histo1D({"dis", "", N_xq2bins, 0, N_xq2bins}, "bin_xBQ2_Valerii");
    s9->Write();
    tfn->Close();
    cout << sname << endl;
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
