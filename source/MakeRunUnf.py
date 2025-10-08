import subprocess
import sys

# Check if the argument is provided
if len(sys.argv) != 2:
    print("Usage: python script.py <i_value> ()")
    sys.exit(1)

# Get # of batch size from the command-line argument
i_batch = sys.argv[1]
banch_size = 10

i = banch_size*i_batch

include_paths = [
    "/u/scigroup/cvmfs/hallb/clas12/sw/noarch/clas12-qadb/3.1.0/srcC/include/",
    "/u/scigroup/cvmfs/hallb/clas12/sw/noarch/clas12-qadb/3.1.0/srcC/rapidjson/include"
]

# Construct the ROOT command script
root_commands = f"""
.I {include_paths[0]}
.I {include_paths[1]}
.L source/make_rec_unfold_matricies.cxx
make_rec_unfold_matricies({i},{banch_size})
"""

# Full shell command
shell_command = f'''
csh
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12root
module load root
root -l -q <<< "{root_commands}"
'''

# Run everything in a shell
subprocess.run(shell_command, shell=True)