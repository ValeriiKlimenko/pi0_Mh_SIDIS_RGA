// g++ -O2 -std=c++17 `root-config --cflags --libs` make_dis_histos_per_file.cxx -o make_dis_histos_per_file
//
// Usage:
//   ./make_dis_histos_per_file <do_data> <do_rec> <do_gen>
//
// where <do_*> is: true/false, 1/0, yes/no, y/n
//
// Hardcoded placeholder directories:
//   - DIS Data  input: data_in_dir
//               output: data_out_dir
//   - DIS Rec   input: rec_in_dir
//               output: rec_out_dir
//   - DIS Gen   input: gen_in_dir
//               output: gen_out_dir
//
// For each input file it creates a *separate* output file containing:
//   Dis_Data:  TH1D h_xQ2_data       (bin_xBQ2_Valerii)
//   Dis_Rec:   TH1D h_xQ2_rec        (bin_xBQ2_Valeriim)
//              TH2D h_xQ2_response   (bin_xBQ2_Valerii vs bin_xBQ2_Valeriim)
//   Gen_Dis:   TH1D h_xQ2_gen        (bin_xBQ2_Valerii)

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;

// -------- configuration (PLACEHOLDER PATHS) -----------------

// Replace these with your actual directories:
const fs::path data_in_dir  = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_dis/dis_data_added_binning";
const fs::path rec_in_dir   = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec_dis/dis_rec_added_binning";
const fs::path gen_in_dir   = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec_dis/dis_gen_added_binning";

const fs::path data_out_dir = "unfolding_dis/";
const fs::path rec_out_dir  = "unfolding_dis/";
const fs::path gen_out_dir  = "unfolding_dis/";

// ------------------------------------------------------------

static const int    nX   = 21;   // xq2bin: 0..20
static const double x_lo = -0.5;
static const double x_hi = x_lo + nX; // 20.5

// Collect all .root files directly in 'dir' (non-recursive)
void collect_root_files(const fs::path &dir, std::vector<fs::path> &out) {
    if (!fs::exists(dir) || !fs::is_directory(dir)) {
        std::cerr << "WARNING: directory does not exist or is not a directory: "
                  << dir << "\n";
        return;
    }

    for (const auto &entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        const auto &p = entry.path();
        if (p.extension() == ".root") {
            out.push_back(p);
        }
    }
}

// Parse a bool from string: true/false, 1/0, yes/no, y/n
bool parse_bool_arg(const std::string &s) {
    std::string l = s;
    std::transform(l.begin(), l.end(), l.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (l == "1" || l == "true"  || l == "yes" || l == "y") return true;
    if (l == "0" || l == "false" || l == "no"  || l == "n") return false;

    std::cerr << "WARNING: could not interpret '" << s
              << "' as boolean, defaulting to false\n";
    return false;
}

// ---------------------- processing functions -------------------------

void process_dis_data_files(const fs::path &in_dir, const fs::path &out_dir) {
    std::vector<fs::path> files;
    collect_root_files(in_dir, files);

    if (files.empty()) {
        std::cout << "No DIS Data files found in " << in_dir << "\n";
        return;
    }

    // Convert paths to strings for RDataFrame
    std::vector<std::string> file_names;
    file_names.reserve(files.size());
    for (const auto &p : files) {
        file_names.push_back(p.string());
    }

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "ERROR: could not create output dir " << out_dir
                  << " : " << ec.message() << "\n";
        return;
    }

    // Single combined output file for all DIS data
    fs::path out_path = out_dir / "dis_data_hists_all.root";

    std::cout << "Processing " << files.size() << " DIS Data files in " << in_dir
              << " -> " << out_path << "\n";

    // RDataFrame over all input files
    ROOT::RDataFrame df("h22", file_names);

    TFile out(out_path.string().c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: cannot create output file " << out_path << "\n";
        return;
    }
    out.cd();

    auto h_xQ2_data = df.Histo1D(
        {"h_xQ2_data",
         "DIS Data (all files combined);bin_xBQ2_Valerii;Events",
         nX, x_lo, x_hi},
        "bin_xBQ2_Valerii"
    );

    h_xQ2_data->Write();
    out.Close();
}

void process_dis_rec_files(const fs::path &in_dir, const fs::path &out_dir) {
    std::vector<fs::path> files;
    collect_root_files(in_dir, files);

    if (files.empty()) {
        std::cout << "No DIS Rec files found in " << in_dir << "\n";
        return;
    }

    std::vector<std::string> file_names;
    file_names.reserve(files.size());
    for (const auto &p : files) {
        file_names.push_back(p.string());
    }

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "ERROR: could not create output dir " << out_dir
                  << " : " << ec.message() << "\n";
        return;
    }

    // Single combined output file for all DIS rec
    fs::path out_path = out_dir / "dis_rec_hists_all.root";

    std::cout << "Processing " << files.size() << " DIS Rec files in " << in_dir
              << " -> " << out_path << "\n";

    ROOT::RDataFrame df("h22", file_names);

    TFile out(out_path.string().c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: cannot create output file " << out_path << "\n";
        return;
    }
    out.cd();

    // 1D reco histogram: bin_xBQ2_Valeriim (all files combined)
    auto h_xQ2_rec = df.Histo1D(
        {"h_xQ2_rec",
         "DIS Rec (measured, all files combined);bin_xBQ2_Valeriim;Events",
         nX, x_lo, x_hi},
        "bin_xBQ2_Valeriim"
    );

    // 2D response matrix: true vs measured bins (all files combined)
    auto h_xQ2_response = df.Histo2D(
        {"h_xQ2_response",
         "DIS Response Matrix (all files combined);true bin_xBQ2_Valerii;reco bin_xBQ2_Valeriim",
         nX, x_lo, x_hi,
         nX, x_lo, x_hi},
        "bin_xBQ2_Valerii",   // x-axis: true
        "bin_xBQ2_Valeriim"   // y-axis: reco
    );

    h_xQ2_rec->Write();
    h_xQ2_response->Write();
    out.Close();
}

void process_dis_gen_files(const fs::path &in_dir, const fs::path &out_dir) {
    std::vector<fs::path> files;
    collect_root_files(in_dir, files);

    if (files.empty()) {
        std::cout << "No DIS Gen files found in " << in_dir << "\n";
        return;
    }

    std::vector<std::string> file_names;
    file_names.reserve(files.size());
    for (const auto &p : files) {
        file_names.push_back(p.string());
    }

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "ERROR: could not create output dir " << out_dir
                  << " : " << ec.message() << "\n";
        return;
    }

    // Single combined output file for all DIS gen
    fs::path out_path = out_dir / "dis_gen_hists_all.root";

    std::cout << "Processing " << files.size() << " DIS Gen files in " << in_dir
              << " -> " << out_path << "\n";

    ROOT::RDataFrame df("h22", file_names);

    TFile out(out_path.string().c_str(), "RECREATE");
    if (out.IsZombie()) {
        std::cerr << "ERROR: cannot create output file " << out_path << "\n";
        return;
    }
    out.cd();

    auto h_xQ2_gen = df.Histo1D(
        {"h_xQ2_gen",
         "DIS Gen (all files combined);bin_xBQ2_Valerii;Events",
         nX, x_lo, x_hi},
        "bin_xBQ2_Valerii"
    );

    h_xQ2_gen->Write();
    out.Close();
}


int run_make_dis_histos(bool do_data, bool do_rec, bool do_gen) {
    std::cout << "Configuration:\n"
              << "  DIS Data: " << (do_data ? "ON" : "OFF") << "\n"
              << "  DIS Rec : " << (do_rec  ? "ON" : "OFF") << "\n"
              << "  DIS Gen : " << (do_gen  ? "ON" : "OFF") << "\n";

    if (do_data) process_dis_data_files(data_in_dir, data_out_dir);
    if (do_rec)  process_dis_rec_files(rec_in_dir,  rec_out_dir);
    if (do_gen)  process_dis_gen_files(gen_in_dir,  gen_out_dir);

    std::cout << "Done.\n";
    return 0;
}

void run_make_dis_histos_root(int do_data = 1, int do_rec = 1, int do_gen = 1) {
    run_make_dis_histos(do_data, do_rec, do_gen);
}
