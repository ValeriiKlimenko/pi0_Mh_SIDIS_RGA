#include <ROOT/RDataFrame.hxx>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TString.h>
#include <vector>
#include <string>
#include <iostream>


#include "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/source/cuts.cxx"
#include "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/source/binning_params.cxx"

using namespace std;

string out_file = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/genStudy/";

// read all files in  a folder
vector<string> getRootFiles(const string& folderPath) {
    vector<string> rootFiles;
    TSystemDirectory dir("dir", folderPath.c_str());
    TList* files = dir.GetListOfFiles();
    if (!files) return rootFiles;

    TIter next(files);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        string fname = file->GetName();
        if (!file->IsDirectory() && fname.find(".root") != string::npos) {
            rootFiles.push_back(folderPath + "/" + fname);
        }
    }
    return rootFiles;
}

// extract subvector for batch analysis
std::vector<std::string> getSubVector(const std::vector<std::string>& input, size_t i, size_t j) {
    size_t end = std::min(i + j, input.size());
    if (i >= input.size()) return {}; // Return empty if i is out of bounds
    return std::vector<std::string>(input.begin() + i, input.begin() + end);
}


ROOT::RDF::RNode GetBinning(ROOT::RDF::RNode node)
{
   return node              
              // All generated particles
              .Define("em_P",                                   "return (TLorentzVector){ex, ey, ez, sqrt(ex*ex+ey*ey+ez*ez + m_e*m_e)};")
              .Define("pi0_Pm",                                 "return (TLorentzVector){px, py, pz, sqrt(px*px+py*py+pz*pz + m_pi*m_pi)};");
    /* 
              .Define("qm",                   q,                {"em_P"})
              .Define("Q2m",                  Q2,               {"qm"})
              .Define("xBm",                  xB,               {"Q2m", "qm"})
              .Define("phi_trentom",          phi_trento,       {"em_P", "qm", "pi0_Pm"})
              .Define("pi0_sidis_PTm",                          "return (pi0_Pm.Vect()).Perp(qm.Vect());")
              .Define("pi0_sidis_PT2m",                         "return pow(pi0_sidis_PTm, 2);")
              .Define("photon_Em",                              "return qm.E();")
               // Marshall had pi0_P instead of pi0_Pm, why?
              .Define("zm",                                     "return pi0_Pm.E() / photon_Em;")//Mimics others
     
              .Define("bin_xBQ2_Valeriim",        bin_xBQ2_Valerii,    {"xBm", "Q2m"})
              .Define("zpt2phit_8x8x9m",       zpt2phit_8x8x9,   {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"});
  */

}


int Get_GenPi0(const int i_batch, const int batch_size) {
    // input files
    string folderPath = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/";
    auto files = getRootFiles(folderPath);
    files = getSubVector(files, i_batch,batch_size);
    
    if (files.empty()) {
        cout << "No ROOT files found." << endl;
        return 1;
    }

    // RDF preparation
    ROOT::RDataFrame df("h22", files); // Replace "treeName" with your actual TTree name
    auto r1 = ROOT::RDF::RNode(df);
    r1 = GetBinning(r1);


    //RDF saving:


    TString rname = Form("%sGen_part_%d.root", out_file.c_str(), i_batch);
    auto tfn = new TFile(rname.Data(), "recreate");
    //r1.Snapshot("h22", rname.Data(), {"bin_xBQ2_Valeriim", "zpt2phit_8x8x9m"});  
    r1.Snapshot("h22", rname.Data(), {"em_P"});  
    return 0;
}

