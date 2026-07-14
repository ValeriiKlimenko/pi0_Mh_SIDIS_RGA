// File: plot_bin_migr_hists.cxx
// Usage from ROOT:
//   root -l -q 'plot_bin_migr_hists.cxx("unfolding","plots")'
//
//  - histDir: directory with h3_*.root produced by define_bin_migr_hists
//  - plotDir: directory where PNGs will be written

#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH3.h>
#include <TH2.h>

#include <iostream>
#include <string>

// Helper: create 2D projections and save as PNG
void SaveProjections(TH3 *h3, const std::string &plotDir)
{
  if (!h3) return;

  std::string baseName = h3->GetName();

  // XY projection
  {
    TH2 *h2 = (TH2*) h3->Project3D("xy"); // x vs y
    if (h2) {
      std::string cname = baseName + "_XY";
      TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 800, 600);
      c->SetLogz();
      h2->SetTitle((baseName + " (XY projection)").c_str());
      h2->Draw("COLZ");
      std::string outName = plotDir + "/" + baseName + "_XY.png";
      c->SaveAs(outName.c_str());
      delete c;
      delete h2;
    }
  }

  // XZ projection
  {
    TH2 *h2 = (TH2*) h3->Project3D("xz"); // x vs z
    if (h2) {
      std::string cname = baseName + "_XZ";
      TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 800, 600);
      c->SetLogz();
      h2->SetTitle((baseName + " (XZ projection)").c_str());
      h2->Draw("COLZ");
      std::string outName = plotDir + "/" + baseName + "_XZ.png";
      c->SaveAs(outName.c_str());
      delete c;
      delete h2;
    }
  }

  // YZ projection
  {
    TH2 *h2 = (TH2*) h3->Project3D("yz"); // y vs z
    if (h2) {
      std::string cname = baseName + "_YZ";
      TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), 800, 600);
      c->SetLogz();
      h2->SetTitle((baseName + " (YZ projection)").c_str());
      h2->Draw("COLZ");
      std::string outName = plotDir + "/" + baseName + "_YZ.png";
      c->SaveAs(outName.c_str());
      delete c;
      delete h2;
    }
  }
}

// Main function
void plot_bin_migr_hists(const std::string &histDir = "unfolding",
                         const std::string &plotDir = "plots")
{
  gROOT->SetBatch(kTRUE);  // no GUI popups

  // Make sure output directory exists
  gSystem->MakeDirectory(plotDir.c_str());

  TSystemDirectory dir(histDir.c_str(), histDir.c_str());
  TList *files = dir.GetListOfFiles();

  if (!files) {
    std::cerr << "No files found in: " << histDir << std::endl;
    return;
  }

  TIter next(files);
  TSystemFile *fSys = nullptr;

  while ((fSys = (TSystemFile*) next())) {
    std::string fname = fSys->GetName();

    // Skip directories and non-ROOT files
    if (fSys->IsDirectory()) continue;
    if (fname.size() < 5 || fname.substr(fname.size()-5) != ".root") continue;
    
    if (fname.find("_fitted") != std::string::npos) continue;

    std::string fullName = histDir + "/" + fname;
    std::cout << "Opening file: " << fullName << std::endl;

    TFile f(fullName.c_str(), "READ");
    if (f.IsZombie()) {
      std::cerr << "  ERROR: cannot open " << fullName << std::endl;
      continue;
    }

    // Loop over keys in the file
    TIter nextKey(f.GetListOfKeys());
    TKey *key = nullptr;

    while ((key = (TKey*) nextKey())) {
      TObject *obj = key->ReadObj();
      if (!obj) continue;

      // Case 1: TH3 stored at top level
      if (obj->InheritsFrom(TH3::Class())) {
        TH3 *h3 = (TH3*) obj;
        std::cout << "  Found TH3: " << h3->GetName() << std::endl;
        SaveProjections(h3, plotDir);
      }

      // Case 2: there is a directory (as in your writer code)
      else if (obj->InheritsFrom(TDirectory::Class())) {
        TDirectory *subdir = (TDirectory*) obj;
        std::cout << "  Entering directory: " << subdir->GetName() << std::endl;

        TIter nextKey2(subdir->GetListOfKeys());
        TKey *key2 = nullptr;
        while ((key2 = (TKey*) nextKey2())) {
          TObject *obj2 = key2->ReadObj();
          if (!obj2) continue;
          if (!obj2->InheritsFrom(TH3::Class())) continue;

          TH3 *h3 = (TH3*) obj2;
          std::cout << "    Found TH3: " << h3->GetName() << std::endl;
          SaveProjections(h3, plotDir);
        }
      }
    }

    f.Close();
  }

  std::cout << "Done. Plots are in: " << plotDir << std::endl;
}
