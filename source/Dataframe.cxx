using namespace std;

#include "/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/source/cuts.cxx"

 // Mars's one:
 //auto poly_mars = makeTH2PolyMars();
double pi_n = TMath::Pi();
//////////////////////////////////////// ANALYSIS: /////////////////////////////////////

// Diffine variables that are the same for data/rec/gen/dis
ROOT::RDF::RNode AddDefine_Kinematics(ROOT::RDF::RNode node)
{

  auto& ctx = get_binning_context();
  
   return node.Define("e_P",                  Get4mom_corr,      {"ex", "ey", "ez", "esec"})
              .Define("e_PT",                                   "return e_P.Pt();")
              .Define("e_mom",                                  "return e_P.P();")
              .Define("e_theta",              theta,            {"e_P"})
              .Define("e_phi",                phi,              {"e_P"})
              //////////////////////////////////////////////////////////////////////
              .Define("q",                    q,                {"e_P"})
              .Define("photon_E",                               "return q.E();")
              .Define("Q2",                   Q2,               {"q"})
              .Define("W2",                   W2,               {"q"})//Mimics others
              .Define("W",                                      "return sqrt(W2);")
              .Define("xB",                   xB,               {"Q2", "q"})
              .Define("y",                    y,                {"q"})//Mimic others
              .Define("gamma",                                  "return 2*m_p*xB / sqrt(Q2);")//Mimics others
              //// Binning: //////////////////
              .Define("bin_xBQ2_Valerii",        ctx.bin_xBQ2,    {"xB", "Q2"});
}
/// Columns for cuts only. They are used in CSB estimation only for now.
ROOT::RDF::RNode AddDefine_CutsCol(ROOT::RDF::RNode node)
{
   return node.Define("strict",            get_cuts_strictness,   {})
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


              .Define("DC_cut",       cut_DC_edge,        {"e_edge_R1","e_edge_R2","e_edge_R3", "strict"}) 
              .Define("SF_cut",       cut_SF,        {"SF_full","e_mom","esec_int", "strict"});

}
////////////// for Both Data and Rec:
ROOT::RDF::RNode AddDefine_Kinematics_RecData(ROOT::RDF::RNode node)
{

  auto& ctx = get_binning_context();

  
   return node.Define("g1_P",                                 "return (TLorentzVector) {g1x, g1y, g1z, sqrt(g1x*g1x+g1y*g1y+g1z*g1z)};")
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
              
              .Define("g1_theta",             theta,          {"g1_P"})
              .Define("g1_phi",               phi,            {"g1_P"})
              .Define("g2_theta",             theta,          {"g2_P"})
              .Define("g2_phi",               phi,            {"g2_P"})
              .Define("pi0_theta",            theta,          {"pi0_P"})
              .Define("pi0_phi",              phi,            {"pi0_P"})

              .Define("g_diff_mom",                            "return (g1_mom-g2_mom);")
              .Define("g_rat_mom",                             "return (g1_mom / g2_mom);")
              .Define("g_prod_mom",                            "return (g1_mom * g2_mom);")
              .Define("g_diff_theta",                          "return (g1_theta-g2_theta);")
              .Define("g_diff_phi",                            "return (g1_phi-g2_phi);")
              .Define("g_open_ang",         ang_bet_vec,       {"g1_P", "g2_P"})
              .Define("e_g1_ang",           ang_bet_vec,       {"e_P", "g1_P"})
              .Define("e_g2_ang",           ang_bet_vec,       {"e_P", "g2_P"})
              .Define("e_pi0_ang",          ang_bet_vec,       {"e_P", "pi0_P"})
              .Define("g_E_balance",                           "return (g1_P.E() - g2_P.E())/(g1_P.E() + g2_P.E());")

              .Define("W2_prime",             W2_prime,         {"q", "pi0_P"})
              .Define("W_prime",                                "return sqrt(W2_prime);")
              .Define("z",                                      "return pi0_P.E() / photon_E;")//Mimics others
              .Define("phi_trento",           phi_trento,       {"e_P", "q", "pi0_P"})
              .Define("Mx2",                  Mx2,              {"q", "pi0_P"})
              .Define("Mx",                                     "return sqrt(Mx2);")
              .Define("ele_ang_bet_vec",      ang_bet_vec,      {"e_P", "q"})
              .Define("pi0_rot",              pi0_rotToQ,       {"q","pi0_P"})
              .Define("pi0_sidis_PT",                           "return (pi0_P.Vect()).Perp(q.Vect());")
              .Define("pi0_sidis_PT2",                          "return pow(pi0_sidis_PT, 2);")
              .Define("xF",                   xF,               {"pi0_P", "q", "W2"})//Mimics the definition on the github as of 10/27/2021
  


     

               // Binning
              // I do not think that update is required it return 10*pt bin + z bin. 2 bins are probably overflow. 
              //.Define("zpt2_8x8",               zpt2_8x8,               {"z", "pi0_sidis_PT2"})
              // I am not sure if it is used anywhere so I will keep it at is for now (dec 18 2024)
              //.Define("xq2zpt2_13x8x8",         xq2zpt2_13x8x8,         {"xB", "Q2", "z", "pi0_sidis_PT2"})
              //Valerii: the binning that was actually used. It is a linearization of multidimensional binning
              // Later on Rebin is used to integrate over phitrento
              .Define("zpt2phit_8x8x9",       ctx.zpt2phit_8x8x9,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})
              .Define("isEventINbins",       ctx.isEventInBins,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})

              // I am not sure if those binngs are used anywhere so I will keep it at is for now (dec 18 2024)
              //.Define("xq2zpt2phit_13x8x8x9",   xq2zpt2phit_13x8x8x9,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})
              //.Define("xq2zpt2phit_13x8x8x12",  xq2zpt2phit_13x8x8x12,  {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"});

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

}

////////////// for Gen Only:
ROOT::RDF::RNode AddDefine_Kinematics_GenOnly(ROOT::RDF::RNode node)
{
  auto& ctx = get_binning_context();
   return node
              .Define("pi0_P",                                "return (TLorentzVector) {px, py, pz, sqrt(px*px+py*py+pz*pz+m_pi*m_pi)};")
              .Define("pi0_PT",                               "return pi0_P.Pt();")
              .Define("pi0_mom",                              "return pi0_P.P();")
              .Define("pi0_m",                                "return pi0_P.M();")
              .Define("pi0_E",                                "return pi0_P.E();")
              .Define("pi0_theta",            theta,          {"pi0_P"})
              .Define("pi0_phi",              phi,            {"pi0_P"})
              .Define("z",                                      "return pi0_P.E() / photon_E;")//Mimics others
              .Define("phi_trento",           phi_trento,       {"e_P", "q", "pi0_P"})
              .Define("Mx2",                  Mx2,              {"q", "pi0_P"})
              .Define("Mx",                                     "return sqrt(Mx2);")
              .Define("ele_ang_bet_vec",      ang_bet_vec,      {"e_P", "q"})
              .Define("pi0_rot",              pi0_rotToQ,       {"q","pi0_P"})
              .Define("pi0_sidis_PT",                           "return (pi0_P.Vect()).Perp(q.Vect());")
              .Define("pi0_sidis_PT2",                          "return pow(pi0_sidis_PT, 2);")
              .Define("xF",                   xF,               {"pi0_P", "q", "W2"})//Mimics the definition on the github as of 10/27/2021
              .Define("isEventINbins",       ctx.isEventInBins,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})
              .Define("zpt2phit_8x8x9",       ctx.zpt2phit_8x8x9,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"});
}
////////////// for Rec Only:

ROOT::RDF::RNode AddDefine_Kinematics_RecOnly(ROOT::RDF::RNode node)
{
  auto& ctx = get_binning_context();
  
   return node              //Matching//////////////////////////////////////////////////////////////
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

              .Define("bin_xBQ2_Valeriim",        ctx.bin_xBQ2,    {"xBm", "Q2m"})

              //.Define("zpt2_8x8m",               zpt2_8x8,               {"zm", "pi0_sidis_PT2m"})
              //.Define("xq2zpt2_13x8x8m",         xq2zpt2_13x8x8,         {"xBm", "Q2m", "zm", "pi0_sidis_PT2m"})
     
               //the one that we need:
              .Define("zpt2phit_8x8x9m",       ctx.zpt2phit_8x8x9,   {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})
              .Define("isEventINbins_m",       ctx.isEventInBins,   {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})
    
              //.Define("xq2zpt2phit_13x8x8x9m",   xq2zpt2phit_13x8x8x9,   {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})
              //.Define("xq2zpt2phit_13x8x8x12m",  xq2zpt2phit_13x8x8x12,  {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})

              .Define("m1_dtheta",                              "return 180.*(g1_P.Theta() - m1_P.Theta())/pi_n;")
              .Define("m2_dtheta",                              "return 180.*(g2_P.Theta() - m2_P.Theta())/pi_n;")
              .Define("m1_dphi",                                "double dphi = g1_P.Phi() - m1_P.Phi(); if (dphi > pi_n){dphi-=2.*pi_n;}  if (dphi < -pi_n){dphi+=2.*pi_n;} return dphi*180.0/pi_n;")
              .Define("m2_dphi",                                "double dphi = g2_P.Phi() - m2_P.Phi(); if (dphi > pi_n){dphi-=2.*pi_n;}  if (dphi < -pi_n){dphi+=2.*pi_n;} return dphi*180.0/pi_n;")
              .Define("m_dtheta_g",                             "return 180.*(m1_P.Theta() - m2_P.Theta())/pi_n;")
              .Define("m_dphi_g",                               "double dphi = m1_P.Phi() - m2_P.Phi(); if (dphi > pi_n){dphi-=2.*pi_n;}  if (dphi < -pi_n){dphi+=2.*pi_n;} return dphi*180.0/pi_n;")
              .Define("m1_dp",                                  "return g1_P.P() - m1_P.P();")
              .Define("m2_dp",                                  "return g2_P.P() - m2_P.P();")
              .Define("m_gmean",                                "return sqrt(g1mVal*g2mVal);")
              .Define("m_gratio",                               "return g1mVal/g2mVal;")
              .Define("nlog_g1mVal",                            "-TMath::Log10(g1mVal);")
              .Define("nlog_g2mVal",                            "-TMath::Log10(g2mVal);")
              .Define("dz",                                     "return z-zm;");
}
////////////// for Unfolding file Only:

/*
ROOT::RDF::RNode AddDefine_Kinematics_UnfOnly(ROOT::RDF::RNode node)
{
   return node//.Define("zpt2_8x8",               zpt2_8x8,               {"z", "pi0_sidis_PT2"})
              //.Define("xq2zpt2_13x8x8",         xq2zpt2_13x8x8,         {"xB", "Q2", "z", "pi0_sidis_PT2"})
              //.Define("zpt2phit_8x8x9",       zpt2phit_8x8x9,   {"z", "pi0_sidis_PT2", "phi_trento"})
              .Define("zpt2phit_8x8x9",       zpt2phit_8x8x9,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})
              .Define("bin_xBQ2_Valerii",        bin_xBQ2_Valerii,    {"xB", "Q2"})

              //.Define("xq2zpt2phit_13x8x8x9",   xq2zpt2phit_13x8x8x9,   {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})
              //.Define("xq2zpt2phit_13x8x8x12",  xq2zpt2phit_13x8x8x12,  {"xB", "Q2", "z", "pi0_sidis_PT2", "phi_trento"})

              //Matching//////////////////////////////////////////////////////////////
              .Define("bin_xBQ2_Valeriim",        bin_xBQ2_Valerii,    {"xBm", "Q2m"})

              //.Define("zpt2_8x8m",               zpt2_8x8,               {"zm", "pi0_sidis_PT2m"})
              //.Define("xq2zpt2_13x8x8m",         xq2zpt2_13x8x8,         {"xBm", "Q2m", "zm", "pi0_sidis_PT2m"})
              //.Define("zpt2phit_8x8x9m",       zpt2phit_8x8x9,   {"zm", "pi0_sidis_PT2m", "phi_trentom"})
              .Define("zpt2phit_8x8x9m",       zpt2phit_8x8x9,   {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})

              //.Define("xq2zpt2phit_13x8x8x9m",   xq2zpt2phit_13x8x8x9,   {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})
    
              //.Define("xq2zpt2phit_13x8x8x12m",  xq2zpt2phit_13x8x8x12,  {"xBm", "Q2m", "zm", "pi0_sidis_PT2m", "phi_trentom"})
              ;
           


}

*/