#!/usr/bin/env python3
import argparse
import os
import re
import shlex
import subprocess
import time
from pathlib import Path

# ----------------- CONFIG: set your prepared data directories here -----------------
DATA_DIR_MAP = {
    # Adjust these three paths as needed:
    "DATA": "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/",
    "REC" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/",
    "GEN" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/",

    "DIS_DATA" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_dis/",
    "DIS_REC" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec_dis/",
    "DIS_GEN" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec_dis/",
  
}
# ----------------------------------------------------------------------------------

MAX_SCREENS = 35
poll_sec = 1
SCREEN_NAME = "addBinning"

def screen_ls():
    try:
        out = subprocess.check_output(["screen", "-ls"], text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        out = e.output

    sessions = {}
    for line in out.splitlines():
        # lines look like: "    12345.job7   (Detached)"
        m = re.search(r"\s*(\d+)\.([^\s]+)\s+\(([^)]+)\)", line)
        if m:
            pid, name, status = int(m.group(1)), m.group(2), m.group(3)
            if SCREEN_NAME in name:
                sessions[name] = {"pid": pid, "status": status}
    return sessions

def start_screen_session(name: str, cmd: str, log_path: Path = None, workdir: Path = None):
    """
    Launch a detached screen session running `cmd`.
    The session will EXIT when `cmd` finishes (because we use `exec`).
    Output is appended to log_path if provided.
    """
    if log_path:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        command = f"exec {cmd} >> {shlex.quote(str(log_path))} 2>&1"
    else:
        command = f"exec {cmd}"

    subprocess.run(
        ["screen", "-dmS", name, "bash", "-lc", command],
        cwd=str(workdir) if workdir else None,
        check=True,
    )

if __name__ == "__main__":
    # --- CLI args ---
    parser = argparse.ArgumentParser(
        description="Launch add_binning ROOT macro in parallel screen sessions."
    )

    parser.add_argument(
        "--data-type",
        choices=["Data", "Rec", "Gen", "DATA", "REC", "GEN", "Dis_data", "Dis_rec", "DIS_GEN"],
        default="Data",
        help="Which dataset to process: Data, Rec, or Gen (case-insensitive).",
    )
    parser.add_argument(
        "--is-true-gen-event",
        type=int,
        choices=[0, 1],
        default=1,
        help="Pass 0 or 1 to the macro as the is_true_gen_event flag (default: 0).",
    )

    parser.add_argument(
        "--mx-cut-mode",
        type=int,
        choices=[0, 1, 2],
        default=2,
        help="Mx cut at No cut, 1.0, and 1.5 respectively.",
    )

    parser.add_argument(
        "--bin-test",
        type=int,
        choices=[0, 1],
        default=0,
        help="fills z, pt2, phi in the output ttree (0/1).",
    )

    parser.add_argument(
        "--isGoodElec",
        type=int,
        choices=[0, 1],
        default=0,
        help="Requires to have detected elec for SIDIS acceptance.",
    )

    parser.add_argument(
        "--isMCParentCut",
        type=int,
        choices=[0, 1],
        default=0,
        help="Gammas has to be originated from pi0, MC only.",
    )

    parser.add_argument(
        "--noPhiBinning",
        type=int,
        choices=[0, 1],
        default=0,
        help="Integrate over phi.",
    )

    parser.add_argument(
        "--logs-dir",
        default="logs_add_binning",
        help="Directory to store screen logs.",
    )
    args = parser.parse_args()

    # Normalize data type and pick directory/tag
    data_tag = args.data_type.upper()  # "DATA", "REC", "GEN"
    if data_tag not in DATA_DIR_MAP:
        raise SystemExit(f"[ERROR] Unknown data type '{args.data_type}'. Expected one of: Data, Rec, Gen.")
    data_dir = DATA_DIR_MAP[data_tag]

    # Validate input dir
    mypath = Path(data_dir)
    if not mypath.exists() or not mypath.is_dir():
        raise SystemExit(f"[ERROR] Data directory does not exist or is not a directory: {mypath}")

    # Gather input files
    rec_files = [p.name for p in mypath.iterdir() if p.is_file()]
    if not rec_files:
        raise SystemExit(f"[ERROR] No input files found in: {mypath}")

    logs_dir = Path(args.logs_dir)
    is_true_gen_event = args.is_true_gen_event
    mx_cut_mode = int(args.mx_cut_mode)
    is_bin_test = int(args.bin_test)
    isGoodElec = int(args.isGoodElec)
    isMCParentCut = int(args.isMCParentCut)
    noPhiBinning = int(args.noPhiBinning)

    LD_EXPORT = (
        'export LD_LIBRARY_PATH='
        '"/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.9.0/4.3.0/lib64:'
        '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.7/lib:'
        '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.36.04/lib:'
        '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/python/3.13.7/lib:'
        '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.3.0/lib:'
        '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib64:'
        '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib:${LD_LIBRARY_PATH-}"'
    )


    for i_file, file in enumerate(rec_files):
        # wait until there is an available screen slot
        while len(screen_ls()) >= MAX_SCREENS:
            time.sleep(poll_sec)

        file_path = str(mypath / file)
        name = f"{SCREEN_NAME}_{i_file}"           # unique screen session name
        log = logs_dir / f"{name}.log"             # capture stdout/stderr

        # Build the macro call (numeric flag unquoted). Pass the selected data tag.
        # Signature: add_binning.cxx("<file>", "<DATA|REC|GEN>", <is_true_gen_event>)

      
        macro_call = f'source/add_binning.cxx("{file_path}", "{data_tag}", {is_true_gen_event}, {mx_cut_mode}, {is_bin_test}, {isGoodElec}, {isMCParentCut}, {noPhiBinning})'
        root_call = f"root -l -b -q {shlex.quote(macro_call)}"

        cmd = f"""
        echo "Start: $(date)"
        module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
        module load clas12root
        module load root
        {LD_EXPORT}
        {root_call}
        echo "Done: $(date)"
        """

        start_screen_session(name, cmd, log_path=log)
        print(f"Launched {name} [data_type={data_tag}] (is_true_gen_event={is_true_gen_event})  (log: {log})")