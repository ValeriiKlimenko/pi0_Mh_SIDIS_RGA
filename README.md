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

./run_rec
./run_data
./run_gen
./run_dis
./run_gen_dis

# 2. Add binning

 Used now (true gen events)
python run_add_binning.py --data-type Rec --is-true-gen-event 1

 Not used at the moment (has to be updated, for fakes)
python run_add_binning.py --is-true-gen-event 0

# 3. Define bin migration

Default: uses rec_true/
python run_define_bin_migr.py

 Not used at the moment
python run_define_bin_migr.py --dir rec_fake/

# 4. Split and fit (SIM)

 SIM (rec_true)
python run_split_and_fit_unif.py --subdir unfolding_rec_true

# Outdated:

rm -rf /lustre24/expphy/volatile/clas12/valerii/multi_pi0/pi0_mass_fits/unfolding_rec_*
python run_split_and_fit.py


#Not used at the moment

python run_split_and_fit.py --subdir unfolding_rec_fake


# Data and Gen
# 1. Add binning and apply the cuts:

  python run_add_binning.py --data-type Data
  python run_add_binning.py --data-type Gen

## 2-3. Fill Rec Data Hist and fit pi0

  python run_define_bin_migr_data.py

# new step: turn gen ttree into th3f xQ2_gen, z_pt2_phi_gen, nPions_Gen.
  make_gen_binning_2D.cxx

# DATA
  python run_split_and_fit_unif.py --subdir unfolding_rec_data

# Outdated:  
  python run_split_and_fit_data.py

## 4. construct response object and data vector + unfold it:

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