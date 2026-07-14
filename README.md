1) Fit in P theta bins
2) update sigma (3)











# TO do:

1) switch to clas12root GEN/REC analysis

# CLAS12 π⁰ Unfolding — Analysis README
---------

  # 2026 (No GregAI):
  kill the screens:
  screen -ls | grep Detached | awk '{print $1}' | xargs -I {} screen -X -S {} quit

  # 0. Environment setup 
  module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
  module load clas12root
  module load root
  module load iguana
  setenv LD_LIBRARY_PATH /home/valerii/work_disk/RooUnfold-master/build:$LD_LIBRARY_PATH

  # 1. Convert HIPO → ROOT
  python run_rec.py --is_greg_ai 0
  python run_gen.py --base-index 0
  python run_gen.py --base-index 1
  ./run_data.sh --is_greg_ai 0
  
  # 2. Add binning
  python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1 --gg_mom_cut 0
  python run_add_binning.py --data-type Data --mx-cut-mode 2 --gg_mom_cut 0
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2

  # for testing
  python run_add_binning.py --data-type Data --mx-cut-mode 2 --gg_mom_cut 0 --custom-output 1
  

  # 3. Define bin migration
  python run_define_bin_migr_data.py --dir rec_data_05 --mx-cut-mode 2
  python run_define_bin_migr.py --mx-cut-mode 2

  make_gen_binning_2D.cxx
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning", "gen_binning_2D_no_gemc512.root")

  # 4. Split and fit 
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_true
  python run_split_and_fit_unif.py --logic data --subdir unfolding_rec_data_05

# DIS:
---------

  # 1. Convert HIPO → ROOT
  # includes both rec-gen subset and full gen subset
  python run_rec_dis.py
  ./run_dis.sh

  # 2. Add binning
  python run_add_binning.py --data-type Dis_data  --mx-cut-mode 0 
  python run_add_binning.py --data-type Dis_rec  --mx-cut-mode 0 
  python run_add_binning.py --data-type DIS_GEN  --mx-cut-mode 0 

  
  get_dis_hist_for_unf.cxx
  .L ~/work_disk/RooUnfold-master/examples/RooUnfoldExample.cxx
  .L source/unfolding_dis_1D.cxx++
  
# Combined Acceptance:
---------
  # 0. - 1. Are the same (No need for rerun)
  # 2. Add binning
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2 --isGoodElec 1

  # 3. Define bin migration
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_sidis_goodElec", "gen_binning_2D_isGoodElec_no_gemc512.root")

  # 0. Updated Gen workflow:
  
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2 --no_cuts_gen 1


# Greg AI:
---------

  # 1. Convert HIPO → ROOT
  (twice for all the files):
  python run_rec.py --is_greg_ai 1
  ./run_data.sh --is_greg_ai 1

  # 2. Add columns, perfrom cuts:
  python run_add_binning.py --data-type Data --mx-cut-mode 2 --gg_mom_cut 1 --gregs_ai_on 1
  python run_add_binning.py --data-type Rec --mx-cut-mode 2 --gg_mom_cut 1 --gregs_ai_on 1
  python run_add_binning.py --data-type Rec --mx-cut-mode 2 --gg_mom_cut 0 --gregs_ai_on 1
   
  # no need to rerun Gen since it have no idea about gammas.

  # 3. Define bin migration
  python run_define_bin_migr_data.py --dir "../data_data_greg_ai/data_photon_pbtp_035"
  python run_define_bin_migr.py  --dir "../data_rec_greg_ai/rec_photon_pbtp_gg_035"
  python run_define_bin_migr.py  --dir "../data_rec_greg_ai/rec_photon_pbtp_gg_05"
  
  # 4. Split and No(!) fit 
  python run_split_and_fit_unif.py --logic data --subdir unfolding_rec_data_035
  python run_split_and_fit_unif.py --logic data --subdir unfolding_rec_data_05

  python run_split_and_fit_unif.py --logic data --subdir unfolding_data_photon_pbtp_05


  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_photon_pbtp_gg_035
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_photon_pbtp_gg_05

  # 5. 

  
  # old:
 ---------
 ---------
 ---------
 ---------

  # 2. Add binning
  # Used now (true gen events) Mx cut at 1.5 DEFAULT:

  # Mx cut at 1
  python run_add_binning.py --data-type Rec  --mx-cut-mode 1 --is-true-gen-event 1

  # Not used at the moment (has to be updated, for fakes)
  python run_add_binning.py --is-true-gen-event 0

  # 3. Define bin migration
  Default: uses rec_true/
  python run_define_bin_migr.py --mx-cut-mode 2

  # Mx cut at 1
  python run_define_bin_migr.py --mx-cut-mode 1
 
  # Not used at the moment
  python run_define_bin_migr.py --dir rec_fake/

  # 4. Split and fit (SIM)
  # SIM (rec_true)
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_true
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_noMxcut

  # Mx cut at 1
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_Mxcut_1

  # Outdated:
  rm -rf /lustre24/expphy/volatile/clas12/valerii/multi_pi0/pi0_mass_fits/unfolding_rec_*
  python run_split_and_fit.py 

  # Not used at the moment
  python run_split_and_fit.py --subdir unfolding_rec_fake

  # Data and Gen
  # 1. Add binning and apply the cuts:

  # Mx cut at 1.5
  python run_add_binning.py --data-type Data --mx-cut-mode 2
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2
  
  # for acceptance study:
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2 --isGoodElec 1
  
  # Marshall's cut:
  python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1 --isMCParentCut 1  
  python run_define_bin_migr.py --mx-cut-mode 2 --isMCParentCut 1
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_match

  # No phi:
  python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1 --noPhiBinning 1  
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2 --noPhiBinning 1  
  python run_add_binning.py --data-type Data  --mx-cut-mode 2 --noPhiBinning 1  
  python run_define_bin_migr_data.py --dir data_no_phi --noPhiBinning 1

  # test with M:
  .L Generate_Data_noPhi.cxx+
  make_meas_data_only_driver("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/unfolding_data_no_phi/h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root")

  # Use data one with the different path because there is no gen info:
  python run_define_bin_migr_data.py --dir rec_no_phi --noPhiBinning 1 --noPhiBinning_forMC 1
  # Data logic is used for both sets because of the same format:
  python run_split_and_fit_unif.py --logic data --subdir unfolding_data_no_phi --noPhiBinning 1
  python run_split_and_fit_unif_Marshall_2026.py --logic data --subdir unfolding_data_no_phi --noPhiBinning 1
  python run_split_and_fit_unif.py --logic data --subdir unfolding_rec_no_phi --noPhiBinning 1
  python run_split_and_fit_unif_Marshall_2026.py --logic data --subdir unfolding_rec_no_phi --noPhiBinning 1
  
  # Acceptance corrected distribution:
  .L source/acc_corrected_noPhi.cxx
  make_data_gen_over_mc("unfolding_data_no_phi/h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root","unfolding_rec_no_phi/h3_bin_xBQ2_Valerii__zpt2phit_8x8x9__pi0_m_fitted.root","gen_binning_noPhi.root")
  # Move file to a folder then
  .L plot_pt2_SIDIS_acc.cxx+
  plot_by_xq2_z_pt2_nophi("data_times_gen_over_mc.root")

  # Mx cut at 1.0
  python run_add_binning.py --data-type Data --mx-cut-mode 1
  python run_add_binning.py --data-type Gen  --mx-cut-mode 1  

  # 2-3. Fill Rec Data Hist and fit pi0
  python run_define_bin_migr_data.py --mx-cut-mode 2

  # Mx cut at 1.0
  python run_define_bin_migr_data.py --mx-cut-mode 1

  # new step: turn gen ttree into th3f xQ2_gen, z_pt2_phi_gen, nPions_Gen.
  make_gen_binning_2D.cxx

  # Mx cut at 1.5:
  # default parametrs are the same:
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning", "gen_binning_2D.root")
  # no Mx cut:
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_noMxcut", "gen_binning_2D_noMxCut.root")
  # Mx cut at 1.0:
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_Mxcut_1", "gen_binning_2D_MxCut_1.root")
  # for acceptance testing:
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_sidis_goodElec", "gen_binning_2D_isGoodElec.root")
  # for phi fit dependence testing:
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning_no_phi", "gen_binning_noPhi.root", 1)
  
  # DATA
  python run_split_and_fit_unif.py --subdir unfolding_rec_data
  python run_split_and_fit_unif.py --logic data --subdir unfolding_data_noMxcut
  python run_split_and_fit_unif.py --logic data --subdir unfolding_data_Mxcut_1

  # Outdated:  
  python run_split_and_fit_data.py

  ## 4. construct response object and data vector + unfold it:

  #one-pass unfold:

  root [0] .L ~/work_disk/RooUnfold-master/examples/RooUnfoldExample.cxx
  root [1] .L source/onepass_unfold_useGenHist.cxx 
  root [2] onepass_unfold("bbb")
  root [2] onepass_unfold("bayes")
  onepass_unfold("bbb_roo")

  # 2026 Jan version:
  root -l >& out.log
  .L ~/work_disk/RooUnfold-master/examples/RooUnfoldExample.cxx
  .L onepass_unfold_useGenHist_1D.cxx
  onepass_unfold("bayes")
  
  ## Dis:
  
  # after root files have been produced, we can add binning and prepare response object
  # mx-cut-mode willbe set to 0 in the code for DIS so any value can be passed
  python run_add_binning.py --data-type Dis_data  --mx-cut-mode 0 
  python run_add_binning.py --data-type Dis_rec  --mx-cut-mode 0 
  python run_add_binning.py --data-type DIS_GEN  --mx-cut-mode 0 
  get_dis_hist_for_unf.cxx
  .L source/unfolding_dis_1D.cxx++
  
  #Greg's AI:
  python run_rec.py --is_greg_ai 1
  python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1 --gregs_ai_on 1 --gg_mom_cut 1
  python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1 --gregs_ai_on 1 --gg_mom_cut 0
  python run_add_binning.py --data-type Data --mx-cut-mode 2 --gregs_ai_on 1 --gg_mom_cut 1
  python run_add_binning.py --data-type Data --mx-cut-mode 2 --gregs_ai_on 1 --gg_mom_cut 0
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2

  python run_define_bin_migr_data.py --dir data_photon_pbtp
  python run_split_and_fit_unif.py --logic data --subdir unfolding_data_photon_pbtp

  python run_define_bin_migr.py --dir rec_photon_pbtp
  python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_photon_pbtp

  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning", "gen_binning_2D_gregAI.root")
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning", "gen_binning_2D_1023files.root.root")
  
  # Outdated:

  run via root (.L):
  create_response_obj.cxx
  - 2D bayes: perform_unfolding.cxx
  - bin-by-bin manual: perform_unfolding_manual_bbb() //("response_out.root","unfold_out_manual_bbb.root", /*include_mc_stat=*/true, /*eps=*/0.0)' 
  - bin-by-bin rooUnfold (diag only) perform_unfolding(false,0, "response_out.root", "unfold_out_bin_by_bin.root");
  fit_phi_unfolded.cxx

  8) fit phi
  
  9) prepare and unfold DIS
  
  10) Construct the Mh and do bin size correction and RC corrections
  
  # Other:
  
  
  ## Quick housekeeping
  
  **Clean up detached GNU Screen sessions**
  
  screen -ls | grep Detached | awk '{print $1}' | xargs -I {} screen -X -S {} quit
  =======
  # Other:
  
  testing:
  
  python run_add_binning.py --data-type Data --bin-test 1
  
  
  python test_binning.py --input_dir /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/rec_true/
  
  ```bash
# Bayes unfolding (nIter=5)  chi2(unfolded vs MC truth) = 5.80922e+14
    ap.add_argument("--a", default="../manual_bbb_isGoodElec/Marshall_vs_MyBBB/dataGenOverRec_method3_sumPhi_all_then_ratio.csv", help="First CSV (m3-like format)")
    ap.add_argument("--b", default="dataGenOverRec_method3_sumPhi_all_then_ratio.csv", help="Second CSV (m3-like format)")