# CLAS12 π⁰ Unfolding — Analysis README

---

# 0. Environment setup 
```bash
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12root
module load root

Add RooUnfold to the environment
setenv LD_LIBRARY_PATH /home/valerii/work_disk/RooUnfold-master/build:$LD_LIBRARY_PATH

# 1. Convert HIPO → ROOT
python run_rec.py
python run_gen.py
./run_data
./run_dis

# includes both rec-gen subset and full gen subset
python run_rec_dis.py

# 2. Add binning
# Used now (true gen events) Mx cut at 1.5
python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1

#Mx cut at 1
python run_add_binning.py --data-type Rec  --mx-cut-mode 1 --is-true-gen-event 1

 Not used at the moment (has to be updated, for fakes)
python run_add_binning.py --is-true-gen-event 0

# 3. Define bin migration

 Default: uses rec_true/
 python run_define_bin_migr.py --mx-cut-mode 2

 #Mx cut at 1
 python run_define_bin_migr.py --mx-cut-mode 1
 

 Not used at the moment
 python run_define_bin_migr.py --dir rec_fake/

# 4. Split and fit (SIM)

 SIM (rec_true)
 python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_true
 python run_split_and_fit_unif.py --logic sim --subdir unfolding_rec_noMxcut

 #Mx cut at 1
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

  # no phi:
  python run_add_binning.py --data-type Rec  --mx-cut-mode 2 --is-true-gen-event 1 --noPhiBinning 1  
  python run_add_binning.py --data-type Gen  --mx-cut-mode 2 --noPhiBinning 1  
  python run_add_binning.py --data-type Data  --mx-cut-mode 2 --noPhiBinning 1  
  
  

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
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_binning","gen_binning_2D.root")
  # no Mx cut
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_noMxcut","gen_binning_2D_noMxCut.root")
  # Mx cut at 1.0
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_Mxcut_1","gen_binning_2D_MxCut_1.root")
  # for testing:
  make_gen_binning_2D("/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/gen_sidis_goodElec","gen_binning_2D_isGoodElec.root")
  

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


## Dis:

# after root files have been produced, we can add binning and prepare response object
# mx-cut-mode willbe set to 0 in the code for DIS so any value can be passed
python run_add_binning.py --data-type Dis_data  --mx-cut-mode 0 
python run_add_binning.py --data-type Dis_rec  --mx-cut-mode 0 
python run_add_binning.py --data-type DIS_GEN  --mx-cut-mode 0 
get_dis_hist_for_unf.cxx
.L source/unfolding_dis_1D.cxx++


#outdated

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