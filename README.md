screen -ls | grep Detached | awk '{print $1}' | xargs -I {} screen -X -S {} quit

jcache get /mss/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/clasdis_rga_* -e vklimenko@anl.gov

========================================================================================
===================ANALYSIS CHAIN:======================================================
========================================================================================

0) Varaiable set up:
1) 
  module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
  module load clas12root
  module load root

  setenv LD_LIBRARY_PATH /home/valerii/work_disk/RooUnfold-master/build:$LD_LIBRARY_PATH

2) HIPO to ROOT: 
  ./run_rec, ./run_data, ./run_gen, ./run_dis, ./run_gen_dis

3) Add binning (twice because one generates HISTS with cut and other without the true matching cut). One is used to call Fill() other is Fake() in response object:

  python run_add_binning.py --data-type Rec --is-true-gen-event 1

  Not used at the moment (has to be updated):
  python run_add_binning.py --is-true-gen-event 0

3) Define bin migration twice again (by defualt rec_true/ will be used), do not run simult. becuase of out of RAM crash (limit may change in the future and it may work simul.):

  
  python run_define_bin_migr.py 

  not used at the moment:
  python run_define_bin_migr.py --dir rec_fake/

4) split and fit, twice again:

  rm -rf  /lustre24/expphy/volatile/clas12/valerii/multi_pi0/pi0_mass_fits/unfolding_rec_*
  python run_split_and_fit.py
  
  not used at the moment:
  python run_split_and_fit.py --subdir unfolding_rec_fake

========================================================================================
=================== Data and Gen :======================================================
========================================================================================

2) Add binning and apply the cuts:

  python run_add_binning.py --data-type Data
  python run_add_binning.py --data-type Gen

3) Fill Rec Data Hist and fit pi0

  python run_define_bin_migr_data.py
  python run_split_and_fit_data.py

5) construct response object and data vector + unfold it:

  run via root (.L):
  create_response_obj.cxx
  perform_unfolding.cxx
  fit_phi_unfolded.cxx


8) fit phi

9) prepare and unfold DIS

10) Construct the Mh and do bin size correction and RC corrections


========================================================================================
=================== other :======================================================
========================================================================================

testing:
python test_binning.py --input_dir /lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/rec_true/