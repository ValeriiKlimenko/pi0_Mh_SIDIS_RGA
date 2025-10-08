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

#include "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/source/binning_params.cxx"
#include "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/source/functions_fitting.cxx"
#include "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/source/Dataframe.cxx"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;
using namespace std;
using namespace RooFit;

const string rootprefix="unf_f2018_";

vector<string> loadMCFiles(const string& filePath, int startIndex, int batchSize, const string& nSidis_path, int maxFiles = 1) {
    vector<string> mcRecFile;
    string prefixRec = filePath + "rec_f2018_";

    ifstream infile(nSidis_path);
    if (!infile) {
        cerr << "Error opening input file" << endl;
        return mcRecFile;
    }

    string line;
    while (getline(infile, line)) {
        size_t pos = line.find("/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/");
        if (pos != string::npos) {
            line.erase(pos, 67);  // Length of the path to erase
        }
        if (line.size() >= 5 && line.substr(line.size() - 5) == ".hipo") {
            line.erase(line.size() - 5);  // Remove ".hipo"
        }

        mcRecFile.emplace_back(prefixRec + line + ".root");
    }

    int endIndex = min(startIndex + batchSize, min(static_cast<int>(mcRecFile.size()), maxFiles));
    return vector<string>(mcRecFile.begin() + startIndex, mcRecFile.begin() + endIndex);
}

void make_rec_unfold_matricies(int stPointMC, int total_files, const string& nSidis_path, const string&  filePath_OUT_makeUnf)
{
  // for SF cuts parametrs isMC set to use MC params  
   bool isMC = true;

  
  // Error Ignore Level Set
  //gErrorIgnoreLevel = kFatal;
  TH1::SetDefaultSumw2();
  ROOT::EnableImplicitMT();
  
  // Setting for graphs
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  //gStyle->SetPadTopMargin(0.05);//14
  gStyle->SetPadLeftMargin(0.1);//14
  gStyle->SetPadRightMargin(0.14);//.1

  string filePath = filePath_IN;  // Set your path
  auto selectedFiles = loadMCFiles(filePath, stPointMC, 1., nSidis_path, total_files);
  ROOT::RDataFrame r0("h22", selectedFiles);


 //Functions ////////////////////////////////////////////////////////////////////////////////
  
 //Generates the xQ2 binning
 auto poly = makeTH2Poly();
 auto poly_rich = makeTH2PolyRich();
 auto poly_mars = makeTH2PolyMars();

 // Generate outlines for the xQ2 binning 
 auto out_mars = makeMarsOutline();
 auto out_rich = makeRichOutline();

  // Returns bin number in the xB vs. Q^2 binning developed by Stefan
  // see Table 4.2 in clas12 Sidis analysis note
  auto bin_xBQ2=[&](const double &x, const double &Q2){
    return poly->FindBin(x, Q2);
 };

  // Returns bin number in the xB vs. Q^2 binning developed by Richard
  // see Table 4.2 in clas12 Sidis analysis note
  auto bin_xBQ2_Rich=[&](const double &x, const double &Q2){
    return poly_rich->FindBin(x, Q2);
 };

  // Returns bin number in the xB vs. Q^2 binning developed by Marshall
  // see Table 4.2 in clas12 Sidis analysis note
  auto bin_xBQ2_Mars=[&](const double &x, const double &Q2){
    return poly_mars->FindBin(x, Q2);
 };

  // Returns magntiude of the momentum
  auto mom=[](vector<double> &p){
    return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  };

  // Returns invariant mass of pi0
  // m =sqrt( (g1 + g2)^2)
  auto mass=[&](vector<double> &g1, vector<double> &g2){
    TLorentzVector g1Vec(g1[0], g1[1], g1[2], g1[3]);
    TLorentzVector g2Vec(g2[0], g2[1], g2[2], g2[3]);
    return (g1Vec + g2Vec).M();
  };

  // Returns q
  // q = e - e'
  auto q=[&](TLorentzVector &e){
    TLorentzVector ei(0, 0, sqrt(b_E*b_E-m_e*m_e), b_E);
    return ei - e;
  };

  // Returns Theta
  auto theta=[](TLorentzVector &p){
    return 180.0 * p.Theta() / TMath::Pi();
  };
 
 // Returns Theta
  auto ang_bet_vec=[](TLorentzVector &p, TLorentzVector &q){
    TVector3 pVec(p[0], p[1], p[2]);
    TVector3 qVec(q[0], q[1], q[2]);
    double theta = pVec.Unit() * qVec.Unit();
    return 180.0 * TMath::ACos(theta) / TMath::Pi();
  };

 // Returns Theta
  auto ang_bet_e_g=[](TLorentzVector &e, TLorentzVector &g1, TLorentzVector &g2){
    TVector3 eVec(e[0], e[1], e[2]);
    TVector3 g1Vec(g1[0], g1[1], g1[2]);
    TVector3 g2Vec(g2[0], g2[1], g2[2]);
    vector <double> theta = {180.0 * TMath::ACos(eVec.Unit() * g1Vec.Unit()) / TMath::Pi(),
                             180.0 * TMath::ACos(eVec.Unit() * g2Vec.Unit()) / TMath::Pi()};
    return theta;
  };

  // Returns Phi
  auto phi=[](TLorentzVector &p){
    double phi = 180.0 * p.Phi() / TMath::Pi();
    if (phi < 0) {phi += 360;}
    return phi;
  };

  //Returns electron sector based upon phi
  auto e_sector=[](const double &phi){
    int sec = 0;
    if (phi < 45 || phi >= 345){sec = 1;}
    if (phi >= 45  && phi < 105){sec = 2;}
    if (phi >= 105 && phi < 165){sec = 3;}
    if (phi >= 165 && phi < 225){sec = 4;}
    if (phi >= 225 && phi < 285){sec = 5;}
    if (phi >= 285 && phi < 345){sec = 6;}
    return sec;
  };

//Returns pi sector based upon phi
  auto p_sector=[](const double &phi){
    int sec = 0;
    if (phi < 15 || phi >= 315){sec = 1;}
    if (phi >= 15  && phi < 75){sec = 2;}
    if (phi >= 75  && phi < 135){sec = 3;}
    if (phi >= 135 && phi < 195){sec = 4;}
    if (phi >= 195 && phi < 255){sec = 5;}
    if (phi >= 255 && phi < 315){sec = 6;}
    return sec;
  };

  // Returns Q^2
  auto Q2=[](TLorentzVector &q){
    return -q.M2();
  };

 // Returns the boost vector that boosts to the photon-proton cm frame
  auto boostVec=[&](TLorentzVector &q){     
    TLorentzVector proton(0, 0, 0, m_p);
    TLorentzVector cms = proton + q;
    return -cms.BoostVector();
  };

  // Rotates electron vector to q frame
  auto e_rotToQ=[&](TLorentzVector &q, TLorentzVector &e){
    TVector3 p(e[0], e[1], e[2]);
    p.RotateUz(q.Vect().Unit());

    TLorentzVector res; 
    res.SetXYZM(p[0], p[1], p[2], m_e);
    return res;
 };

  // Rotates beam vector to q frame
  auto beam_rotToQ=[&](TLorentzVector &q){
    TVector3 p(0, 0, sqrt(b_E*b_E-m_e*m_e));
    p.RotateUz(q.Vect().Unit());
    
    TLorentzVector res; 
    res.SetXYZM(p[0], p[1], p[2], m_e);
    return res;
 };

  // Rotates pion vector to q frame
  auto pi0_rotToQ=[&](TLorentzVector &q, TLorentzVector &pi){
    TVector3 p(pi[0], pi[1], pi[2]);
    p.RotateUz(q.Vect().Unit());

    TLorentzVector res; 
    res.SetXYZM(p[0], p[1], p[2], pi.M());
    return res;
 };

  // Boosts the momentum to the photon-proton cm frame
  // When performed after the rotation functions, the resulting vector is in the sidis frame
  auto sidis_P=[](TVector3 &b, TLorentzVector &p){
    TLorentzVector p_rot(p[0], p[1], p[2], p[3]);
    return p_rot.Boost(b);
  };

  // Same as above just for proton
  auto sidis_P_proton=[&](TVector3 &b){
    TLorentzVector p_rot(0, 0, 0, m_p);
    return p_rot.Boost(b);
  };

  // Returns |Momentum| transverse to the q direction
  auto pT_q=[](TLorentzVector &p, TLorentzVector &q){
    TVector3 pVec(p[0], p[1], p[2]);
    TVector3 qVec(q[0], q[1], q[2]);
    return (pVec.Cross(qVec.Unit())).Mag();
  };

 // Returns W^2
 // W^2 = (P + q)^2
 // Note that P is proton 4-vector in the lab frame
 auto W2=[&](TLorentzVector &q){
   TLorentzVector proVec(0, 0, 0, m_p);
   return (q + proVec).Mag2();
 };

 // Returns W'^2
 // W'^2 = (m_p + nu - Eh)^2 - (q - Ph)^2
 auto W2_prime=[&](TLorentzVector &q, TLorentzVector &p){
   return pow(m_p + q.E() - p.E(), 2.0) - (q - p).Mag2();
 };

 // Returns W^2 sidis
 // W^2 = (P + q)^2
 auto W2_sidis=[](TLorentzVector &q, TLorentzVector &p){
   return (q + p).Mag2();
 };

 //Returns Bjorken x
 // xB = Q^2 / (2 * P * q)
 // Note that that this is defined in the lab frame, so P is (0, 0, 0, m_p)
 auto xB=[&](const double &Q2, TLorentzVector &q){
  TLorentzVector proVec(0, 0, 0, m_p);
  double num = proVec*q;
  return Q2 / (2.0*num);
 };

 // Returns xF in sidis frame
 // xF = 2*(Ph_vec * q_vec) / (|q_vec| * W)
 // _vec denotes 3 vector
 auto xF=[](TLorentzVector &pi, TLorentzVector &q, const double &W2){
   auto cms = q + (TLorentzVector){0,0,0, m_p};
   auto q_b  = q;  q_b.Boost(-cms.BoostVector());
   auto pi_b = pi; pi_b.Boost(-cms.BoostVector());

   //TVector3 qVec(q[0], q[1], q[2]);
   //TVector3 piVec(pi[0], pi[1], pi[2]);
   double den = pi_b.Vect() * (q_b.Vect()).Unit();
   return 2.0 * den / sqrt(W2);
   //return (pi*q) / (q.Mag() * sqrt(W2));
 };

// Returns y
// y = P*q / P*b , with P being proton, and b being beam lepton 
// All are in the lab frame
 auto y=[&](TLorentzVector &q){
  TLorentzVector pVec(0, 0, 0, m_p);
  TLorentzVector beamVec(0, 0, sqrt(b_E*b_E - m_e*m_e), b_E);
  double den = pVec*q;
  double num = pVec*beamVec;
  return (den / num);
 };

// Returns y in sidis frame
// y = P*q / P*b , with P being proton, and b being beam lepton
// P, b are in sidis frame
 auto y_sidis=[&](TLorentzVector &p, TLorentzVector &q, TLorentzVector &e){
  double den = p*q;
  double num = p*e;
  return (den / num);
 };

// Returns z in sidis frame
// z = P*Ph / P*q, with P being proton, and Ph being the hadron
// P, Ph in sidis frame
  auto z_sidis=[](TLorentzVector &p, TLorentzVector &q, TLorentzVector &pi){
    double den = p*pi;
    double num = p*q;
    return (den / num);
  };

//Returns trento convention phi
// (qxl)*Ph/|(qxl)*Ph| * arccos[((qxl)/|qxl|) * ((qxPh)/|qxPh|)]
// With l being the electron, q the photon, and Ph the hadron
// Defined in the lab frame
 auto phi_trento_check=[](TLorentzVector &e, TLorentzVector &q, TLorentzVector &pi){
   TVector3 eVec(e[0], e[1], e[2]);
   TVector3 qVec(q[0], q[1], q[2]);
   TVector3 piVec(pi[0], pi[1], pi[2]);
   TVector3 piT(pi[0], pi[1] - q[1], pi[2] - q[2]);
   TVector3 vT = qVec.Cross(eVec);
   TVector3 vTH = qVec.Cross(piT); 
   double cosphi = vT.Unit() * vTH.Unit();
   double sinphi = (eVec.Cross(piT)) * qVec;
   double phi = TMath::ACos(cosphi);
   if (sinphi < 0){phi = 2*TMath::Pi() - phi;}
   double deg = 180.0*phi/TMath::Pi();
   return deg;
 };

//Returns trento convention phi
// (qxl)*Ph/|(qxl)*Ph| * arccos[((qxl)/|qxl|) * ((qxPh)/|qxPh|)]
// With l being the electron, q the photon, and Ph the hadron
// Defined in the lab frame
auto phi_trento=[](TLorentzVector &e, TLorentzVector &q, TLorentzVector &pi){
   TVector3 eVec(e[0], e[1], e[2]);
   TVector3 qVec(q[0], q[1], q[2]);
   TVector3 piVec(pi[0], pi[1], pi[2]);
   double cosphi = (qVec.Cross(eVec)).Unit() * (qVec.Cross(piVec)).Unit();
   double sign = (qVec.Cross(eVec)) * piVec;
   sign /= fabs(sign);
   double phi = sign * TMath::ACos(cosphi);
   double deg = 180.0*phi/TMath::Pi();
   if (deg< 0){deg += 360;}
   return deg;
 };

// Returns the missing mass^2
// Mx^2 = (q + P -Ph)^2
// With P being the proton, Ph being the hadron
// All in lab frame
 auto Mx2=[&](TLorentzVector &q, TLorentzVector &pi){
   TLorentzVector proVec(0, 0, 0, m_p);
   return  (q + proVec - pi).Mag2();
};

//Returns true if bars are live and false if dead
 auto deadbar=[&](double &sec, double &lv, double &lw){
     if (sec == 2 && lv > 103.5 && lv < 112.5){return 0;}
     if (sec == 1 && lw > 76.5  && lw < 81){return 0;}
     if (sec == 1 && lw > 85.5  && lw < 90){return 0;}
     return 1;
};

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

  //Recon MonteCarlo DataFrame
  //recrec
  //cout << "right b4 r0 -> r1" << endl;
  auto r1 = r0
              //Momentum, Phi, Theta, Sector : Lab frame /////////////////////////////
              .Define("e_P",                                  "return (TLorentzVector) {ex, ey, ez, sqrt(ex*ex+ey*ey+ez*ez+m_e*m_e)};")
              .Define("e_PT",                                 "return e_P.Pt();")
              .Define("e_mom",                                "return e_P.P();")
              .Define("g1_P",                                 "return (TLorentzVector) {g1x, g1y, g1z, sqrt(g1x*g1x+g1y*g1y+g1z*g1z)};")
              .Define("g1_PT",                                "return g1_P.Pt();")
              .Define("g1_mom",                               "return g1_P.P();")
              .Define("g2_P",                                 "return (TLorentzVector) {g2x, g2y, g2z, sqrt(g2x*g2x+g2y*g2y+g2z*g2z)};")
              .Define("g2_PT",                                "return g2_P.Pt();")
              .Define("g2_mom",                               "return g2_P.P();")
              .Define("pi0_P",                                "return g1_P + g2_P;")
              .Define("pi0_PT",                               "return pi0_P.Pt();")
              .Define("pi0_mom",                              "return pi0_P.P();")
              .Define("pi0_m",                                "return pi0_P.M();")
              .Define("pi0_E",                                "return pi0_P.E();")

              .Define("e_theta",              theta,          {"e_P"})
              .Define("e_phi",                phi,            {"e_P"})
              //.Define("e_sector",           e_sector,         {"e_phi"})
              //.Define("deadbar_ele",          deadbar,        {"esec", "ev", "ew"})


              .Define("g1_theta",             theta,          {"g1_P"})
              .Define("g1_phi",               phi,            {"g1_P"})
              //.Define("g1pcal",                               "return g1pcal;")
              //.Define("g1_sector",          e_sector,         {"e_phi"})
              //.Define("deadbar_gamma1",       deadbar,        {"g1pcal", "g1v", "g1w"}) 

              .Define("g2_theta",             theta,          {"g2_P"})
              .Define("g2_phi",               phi,            {"g2_P"})
              //.Define("g2pcal",                               "return g2pcal;")
              //.Define("g2_sector",          e_sector,         {"e_phi"})
              //.Define("deadbar_gamma2",       deadbar,        {"g2pcal", "g2v", "g2w"}) 

              .Define("pi0_theta",            theta,          {"pi0_P"})
              .Define("pi0_phi",              phi,            {"pi0_P"})
              //.Define("pi0_sector",          e_sector,         {"e_phi"})


              //Differences/////////////////////////////////////////////////////////
              .Define("g_diff_mom",                            "return (g1_mom-g2_mom);")
              .Define("g_rat_mom",                             "return (g1_mom / g2_mom);")
              .Define("g_prod_mom",                            "return (g1_mom * g2_mom);")
              .Define("g_diff_theta",                          "return (g1_theta-g2_theta);")
              .Define("g_diff_phi",                            "return (g1_phi-g2_phi);")
              .Define("g_open_ang",         ang_bet_vec,       {"g1_P", "g2_P"})
              .Define("sin2_half_g_open_ang",                  "return pow(TMath::Sin(g_open_ang*TMath::Pi()/360.0), 2);")//sin^2(ang/2)
              .Define("e_g1_ang",           ang_bet_vec,       {"e_P", "g1_P"})
              .Define("e_g2_ang",           ang_bet_vec,       {"e_P", "g2_P"})
              .Define("e_pi0_ang",          ang_bet_vec,       {"e_P", "pi0_P"})
              .Define("g_E_balance",                           "return (g1_P.E() - g2_P.E())/(g1_P.E() + g2_P.E());")
              .Define("g_diff_v",                              "return (g1v-g2v);")
              .Define("g_diff_w",                              "return (g1w-g2w);")
              //.Define("g_fun",                                 "return func(g_open_ang);")
              //////////////////////////////////////////////////////////////////////
              .Define("q",                    q,                {"e_P"})
              .Define("photon_E",                               "return q.E();")
              .Define("Q2",                   Q2,               {"q"})
              .Define("W2",                   W2,               {"q"})//Mimics others
              .Define("W2_prime",             W2_prime,         {"q", "pi0_P"})
              .Define("W",                                      "return sqrt(W2);")
              .Define("W_prime",                                "return sqrt(W2_prime);")
              .Define("xB",                   xB,               {"Q2", "q"})
              .Define("y",                    y,                {"q"})//Mimic others
              .Define("z",                                      "return pi0_P.E() / photon_E;")//Mimics others
              .Define("gamma",                                  "return 2*0.93827208816*xB / sqrt(Q2);")//Mimics others
              .Define("phi_trento",           phi_trento,       {"e_P", "q", "pi0_P"})
              .Define("Mx2",                  Mx2,              {"q", "pi0_P"})
              .Define("Mx",                                     "return sqrt(Mx2);")
              .Define("ele_ang_bet_vec",      ang_bet_vec,      {"e_P", "q"})
              .Define("pi0_rot",              pi0_rotToQ,       {"q","pi0_P"})
              .Define("pi0_sidis_PT",                           "return (pi0_P.Vect()).Perp(q.Vect());")
              .Define("pi0_sidis_PT2",                          "return pow(pi0_sidis_PT, 2);")
              .Define("xF",                   xF,               {"pi0_P", "q", "W2"})//Mimics the definition on the github as of 10/27/2021 

// bining:
  
              //.Define("bin_xBQ2",             bin_xBQ2,         {"xB", "Q2"})
              //.Define("bin_xBQ2_Rich",        bin_xBQ2_Rich,     {"xB", "Q2"})
              //.Define("bin_xBQ2_Mars",        bin_xBQ2_Mars,    {"xB", "Q2"})

              .Define("bin_xBQ2_Valerii",        bin_xBQ2_Valerii,    {"xB", "Q2"})
  
              //Matching//////////////////////////////////////////////////////////////
              .Define("em_P",                                   "return (TLorentzVector){emx, emy, emz, sqrt(emx*emx+emy*emy+emz*emz + m_e*m_e)};")
              .Define("m1_P",                                   "return (TLorentzVector){g1mx, g1my, g1mz, sqrt(g1mx*g1mx+g1my*g1my+g1mz*g1mz)};")
              .Define("m2_P",                                   "return (TLorentzVector){g2mx, g2my, g2mz, sqrt(g2mx*g2mx+g2my*g2my+g2mz*g2mz)};")
              .Define("pi0_Pm",                                 "return m1_P + m2_P;")
              .Define("pi0_PTm",                                "return pi0_Pm.Pt();")
              .Define("pi0_momm",                               "return pi0_Pm.P();")
              .Define("pi0_mm",                                 "return pi0_Pm.M();")
              .Define("qm",                   q,                {"em_P"})
              .Define("photon_Em",                              "return qm.E();")
              .Define("Q2m",                  Q2,               {"qm"})
              .Define("W2m",                  W2,               {"qm"})//Mimics others
              .Define("W2_primem",            W2_prime,         {"qm", "pi0_Pm"})
              .Define("Wm",                                     "return sqrt(W2m);")
              .Define("W_primem",                               "return sqrt(W2_primem);")
              .Define("xBm",                  xB,               {"Q2m", "qm"})
              .Define("ym",                   y,                {"qm"})//Mimic others
              .Define("zm",                                     "return pi0_P.E() / photon_Em;")//Mimics others
              .Define("gammam",                                 "return 2*0.93827208816*xB / sqrt(Q2m);")//Mimics others
              .Define("phi_trentom",          phi_trento,       {"em_P", "qm", "pi0_Pm"})
              .Define("Mx2m",                 Mx2,              {"qm", "pi0_Pm"})
              .Define("Mxm",                                     "return sqrt(Mx2m);")
              .Define("ele_ang_bet_vecm",     ang_bet_vec,      {"em_P", "qm"})
              .Define("pi0_rotm",             pi0_rotToQ,       {"qm","pi0_Pm"})
              .Define("pi0_sidis_PTm",                          "return (pi0_Pm.Vect()).Perp(qm.Vect());")
              .Define("pi0_sidis_PT2m",                         "return pow(pi0_sidis_PTm, 2);")
              .Define("xFm",                  xF,               {"pi0_Pm", "qm", "W2m"})//Mimics the definition on the github as of 10/27/2021 

// binning update:
  
              //.Define("bin_xBQ2m",            bin_xBQ2,         {"xBm", "Q2m"})
              //.Define("bin_xBQ2_Richm",       bin_xBQ2_Rich,    {"xBm", "Q2m"})
              //.Define("bin_xBQ2_Marsm",       bin_xBQ2_Mars,    {"xBm", "Q2m"})
              .Define("bin_xBQ2_Valeriim",        bin_xBQ2_Valerii,    {"xBm", "Q2m"})

  
              //Valerii params for cuts:
              .Define("strict",            get_cuts_strictness,   {})
              .Define("SF_full",                            "return (e_pcalE + e_ecinE + e_ecoutE)/e_P.P();")
              .Define("SF_pcal_ecin",                            "return (e_pcalE + e_ecinE)/e_P.P();")
              .Define("esec_int",                            "return (int)lrint(esec);")
              .Define("g1sec_int",                            "return (int)lrint(g1sec);")
              .Define("g2sec_int",                            "return (int)lrint(g2sec);")

              //Cuts:
              .Define("cut_pcal_fid_el",       cut_PCAL_fid,        {"e_pcal_Lw", "e_pcal_Lv", "e_pcal_Lu",
                                                                        "e_ecin_Lw", "e_ecin_Lv", "e_ecin_Lu",
                                                                        "e_ecout_Lw", "e_ecout_Lv", "e_ecout_Lu",
                                                                        "esec_int","strict"}) 


              .Define("cut_pcal_fid_g1",       cut_PCAL_fid,        {"g1_pcal_Lw", "g1_pcal_Lv", "g1_pcal_Lu",
                                                                        "g1_ecin_Lw", "g1_ecin_Lv", "g1_ecin_Lu",
                                                                        "g1_ecout_Lw", "g1_ecout_Lv", "g1_ecout_Lu",
                                                                        "g1sec_int","strict"}) 
    
              .Define("cut_pcal_fid_g2",       cut_PCAL_fid,        {"g2_pcal_Lw", "g2_pcal_Lv", "g2_pcal_Lu",
                                                                        "g2_ecin_Lw", "g2_ecin_Lv", "g2_ecin_Lu",
                                                                        "g2_ecout_Lw", "g2_ecout_Lv", "g2_ecout_Lu",
                                                                        "g2sec_int","strict"}) 
              .Define("DC_cut",       cut_DC_edge,        {"e_edge_R1","e_edge_R2","e_edge_R3", "strict"}) 
              .Define("SF_cut",       cut_SF,        {"SF_full","e_mom","esec_int", "strict"});
    
              ;
 

  cout << "right after r0 -> r1" << endl;

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
  //End cuts/////////////////////////////////////////////////////////////////////////


  //Ancillary Cuts
  //Note that for the unfolding matricies, I am dropping any Q2, x, z, pt2 cuts, but only analyzing the events 
  //that make the eventual cut.

  string bin_cuts_dis = "Q2>1 && Q2<12";//Accounting for migration
  //bin_cuts_dis       += "&& xB > 0.1 && xB < 0.8";//Accounting for migration
  string bin_cuts_mul = bin_cuts_dis;

  //bin_cuts_mul       += "&& z > 0.2 && z < 0.8";//0.1-0.8
  //bin_cuts_mul       += "&& pi0_sidis_PT2 > 0 && pi0_sidis_PT2 < 1";
  
  auto r2 = r1.Filter(newCuts.c_str()).Filter(bin_cuts_mul.c_str()).Filter("g1match*g2match>0");
  auto r3 = r2.Filter(tt_cut.c_str()).Filter("pi0_m > 0.1 && pi0_m < 0.168");

  bool compute = 0;
  gROOT->SetBatch(kTRUE);
  
  //Running files
  compute = 1;
  if (compute == 1){

    TString rname = Form("%s%d.root", (filePath_OUT_makeUnf + rootprefix).c_str(), stPointMC);
    rname.ReplaceAll(".hipo", "");
    
    //TString rname = Form("%spi0_rec_unfold_event_%d.root", filePath_OUT_makeUnf.c_str(), start);
    auto tfn = new TFile(rname.Data(), "recreate");
    r2.Snapshot("h22", rname.Data(), {"Mx", "Q2", "xB", "z", "pi0_sidis_PT2", "phi_trento", "pi0_m",
                                      "Mxm","Q2m","xBm","zm","pi0_sidis_PT2m","phi_trentom"});   
  }
}//file end
