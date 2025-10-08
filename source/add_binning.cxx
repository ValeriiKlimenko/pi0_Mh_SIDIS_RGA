// This file reads TTrees into rdf from data, rec and gen root files
// Defines all the variables that are neccessry for the cuts, binning and analysis
// It saves RDF instead with the following columns:  
// For Data 3D: xQ2-bin, z-pt2-phi-Bin, result of 2Y-invariantMass fir
// For Rec vector of histograms:  z-pt2-phi-Bin, 2Y-invariantMass

// g++ -O2 -std=c++17 `root-config --cflags --libs` fit_groups_rdf.cxx -o fit_groups_rdf
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>

#include <unordered_map>
#include <vector>
#include <memory>
#include <tuple>
#include <cmath>
#include <limits>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <string_view>
#include <algorithm>


#include "binning_params.cxx"
#include "Dataframe.cxx"

using ROOT::RDataFrame;
using namespace RooFit;
namespace fs = std::filesystem;









// a flag that save TTree with z information instead of pi0 mass
// the purpose of this to make sure that the bins propogate corretcly
// we have (z,pt2,phi) -> liner binning -> many operations -> unfolding -> unliner to (z,pt2,phi) 
// It is the check that we get the same bin at the very end
bool run_z_bins_instead_of_pi0mass = false;

// this code is for SIDIS part only. There is another file that is used for DIS part.


// Rec means processing for RESPONSE MATRIX, it is not supposed to be used for REC DATA. Use DATA for REC Data.
// I know that it should have been called differentelly, I may fix it later.


enum class FileType {Data, Rec, Gen, Dis, Unknown};

static FileType parse_type_name(string s) {
  auto lower = std::string(s);
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  if (lower == "data") return FileType::Data;
  if (lower == "rec") return FileType::Rec;
  if (lower == "gen") return FileType::Gen;
  if (lower == "dis") return FileType::Dis;
  return FileType::Unknown;
}

void add_binning(const string full_file_path, string TYPE, bool is_true_gen_event) {

  // Type of dara file that is being analyzed
  FileType type;
  type = parse_type_name(TYPE);

  

  // zpPolyVec is globasl so MT causes satbility issues, the number of events with it ON and OFF are very different.
  // it can be fixed since zpPolyVec does not contain actual events and only serves for bin determination. (by making a copy for each slot)
  // however, since it is the only step of the analysis where it is used and it higly parralled (0.25 OSG submission per process)
  // The net weight in the analysis chain is low and not worth the time investment for now.
  
  //ROOT::EnableImplicitMT(); // commented out for binning calculation 

  ////////////////////////////////////////////////////////
  ///// Input Args:
  ////////////////////////////////////////////////////////


    fs::path path_root_in = fs::path(full_file_path);
    fs::path dir = path_root_in.parent_path();

    // set this however you like; using literal "rec" here

    // it was codded before REC data was added so rec is referring to output folder
    fs::path rec;
    if (type == FileType::Rec){
      if (is_true_gen_event) rec = "rec_true";
      else  rec = "rec_fake";
    }
  
    if (type == FileType::Data){
       rec = "rec_data";
    }

    if (type == FileType::Gen){
       rec = "gen_binning";
    }

    // ensure: <dir>/<rec>/<stem>/ exists
    fs::path out_dir = dir / rec;
    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "Failed to create directories: " << out_dir << " : " << ec.message() << "\n";
        return;
    }


  string name_ending = "_out.root";
  if (run_z_bins_instead_of_pi0mass) name_ending = "_out_z_instead_of.root";
  fs::path path_root_out = out_dir / (path_root_in.stem().string() + name_ending);



  if (type != FileType::Rec && type != FileType::Data && type != FileType::Gen) std::cerr << "Wrong file type (not Rec, Data or Gen): " << ' ' << "\n";

  ////////////////////////////////////////////////////////
  ///// Read RDF:
  ////////////////////////////////////////////////////////
  
  bool isMC =  type == FileType::Rec ? true : false;

  // Add the following columns ("xq2bin", "xq2bin_gen", "zpt2bin", "zpt2bin_gen", "MM")
  ROOT::RDataFrame df("h22", path_root_in.c_str());
  auto r1 = ROOT::RDF::RNode(df);

  ////////////////////////////////////////////////////////
  ///// Define columns:
  ////////////////////////////////////////////////////////
  
  r1 = AddDefine_Kinematics(r1);
  // Data and Response Obj:
  if (type != FileType::Gen) r1 = AddDefine_Kinematics_RecData(r1);
  // Response Obj only:
  if (type == FileType::Rec) r1 = AddDefine_Kinematics_RecOnly(r1);
  // Gen only:
  if (type == FileType::Gen) r1 = AddDefine_Kinematics_GenOnly(r1);

  ////////////////////////////////////////////////////////
  ///// Cuts:
  ////////////////////////////////////////////////////////

  // Define all the neccessary columns:
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

  // standard list of cuts Rec, Data, won't work for gen (not enough info)
  // the input parameter means nothing. I may update it to type to use one function for all the cuts.
  string newCuts       = GetMainCuts(true);
  string genCuts       = GetMainCuts_Gen();

  // booking rdf obj to use the same name for muliple datasets
  ROOT::RDF::RNode rdf_after_the_cuts = ROOT::RDF::AsRNode(r1);

  // prepare True Rec and Fake Rec Response rdfs
  if (type == FileType::Rec){
    //MC 4 component cuts
    string sig_cut       = "g1mPID==22 && g2mPID==22 && g1mPPID==111 && g2mPPID==111 && g1mPIndex==g2mPIndex";
    string tt_cut        = sig_cut;
    
    //if (!is_true_gen_event) tt_cut = "!(" + tt_cut + ")";
    //if (!is_true_gen_event) rdf_after_the_cuts = rdf_after_the_cuts.Filter(tt_cut.c_str()); //.Filter("pi0_m > 0.1 && pi0_m < 0.168");
    //rdf_after_the_cuts = rdf_after_the_cuts.Filter(newCuts.c_str()).Filter("g1match*g2match>0 && pi0_sidis_PT2m < 1.5 && isEventINbins_m");
    //rdf_after_the_cuts = rdf_after_the_cuts.Filter(newCuts.c_str()).Filter("g1match*g2match>0 && pi0_sidis_PT2m < 1.5").Filter("isEventINbins");
    
    //rdf_after_the_cuts = rdf_after_the_cuts.Filter(newCuts.c_str()).Filter("isEventINbins").Filter("g1match*g2match>0 && pi0_sidis_PT2m < 1.5");//
    rdf_after_the_cuts = rdf_after_the_cuts.Filter(newCuts.c_str()).Filter("g1match*g2match>0 && pi0_sidis_PT2m < 1.5");//
    
    //rdf_after_the_cuts = rdf_after_the_cuts.Filter("isEventINbins");

    // saved branches depends on dataset, this one is for response matrix:
    rdf_after_the_cuts.Snapshot("h22", path_root_out.string(), {"bin_xBQ2_Valerii","bin_xBQ2_Valeriim","zpt2phit_8x8x9","zpt2phit_8x8x9m","pi0_m","isEventINbins_m","isEventINbins", "xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"});
  }

  // Data Rec distribution. It can be used for filling Rec MC distribution but it is not tested since it is not used in the current workflow. 
  // Mc information comes from response matrix and genrsted rdf.
  if (type == FileType::Data){
    // no new cuts are needed
    rdf_after_the_cuts = rdf_after_the_cuts.Filter(newCuts.c_str());
    // saved branches depends on dataset, this one is for rec data:

    if (!run_z_bins_instead_of_pi0mass) rdf_after_the_cuts.Snapshot("h22", path_root_out.string(), {"bin_xBQ2_Valerii","zpt2phit_8x8x9","pi0_m"});
    if (run_z_bins_instead_of_pi0mass) rdf_after_the_cuts.Snapshot("h22_z", path_root_out.string(), {"bin_xBQ2_Valerii","zpt2phit_8x8x9","z"});
  }
  
  // Gen has unique set of cuts:
  // if an event has multiple gen pi0 it will fill rdf multiple times so there is no need for weight
  if (type == FileType::Gen){
    // no new cuts are needed
    rdf_after_the_cuts = rdf_after_the_cuts.Filter(genCuts.c_str());
    // saved branches depends on dataset, this one is gen:
    // each row is exactly one generated pi0 so there is no need for pi0_m fit.
    // I could have saved hist right away but I will keep intermidaite RDF for consistency with other tyoes of samples.
    rdf_after_the_cuts.Snapshot("h22", path_root_out.string(), {"bin_xBQ2_Valerii","zpt2phit_8x8x9"});
  }

  return;
}
