#!/usr/bin/env python3
import argparse
import os
import re
import shlex
import subprocess
import time
from pathlib import Path

# ----------------- CONFIG -----------------
BASE_PREP_DIR = Path("/lustre24/expphy/volatile/clas12/valerii/multi_pi0")


DATA_DIR_MAP = {
    "GEN":      str(BASE_PREP_DIR / "data_gen"),
    "DIS_DATA": str(BASE_PREP_DIR / "data_dis"),
    "DIS_REC":  str(BASE_PREP_DIR / "data_rec_dis"),
    "DIS_GEN":  str(BASE_PREP_DIR / "data_rec_dis"),
}
# -------------------------------------------

MAX_SCREENS = 30
poll_sec = 1
SCREEN_NAME = "addBinning"


def pick_input_dir(data_tag: str,
                   mx_cut_mode: int,
                   is_true_gen_event: int,
                   is_mc_parent_cuts: int,
                   no_phi_binning: int,
                   gregs_ai_on: int,
                   gg_mom_cut: int,
                   custom_output: int) -> Path:
    """
    Compute input directory name for DATA/REC exactly like the C++ 'rec' label logic.
    Returns an absolute Path.

    For DATA with custom_output enabled, use data_test as input.
    The macro data type remains DATA.
    """
    tag = data_tag.upper()

    # DATA
    if tag == "DATA":
        if custom_output:
            sub = "data_data"
        elif gregs_ai_on:
            sub = "data_data_greg_ai"
        else:
            sub = "data_data"

        return BASE_PREP_DIR / sub

    # REC
    if tag == "REC":
        if gregs_ai_on:
            sub = "data_rec_greg_ai"
        else:
            sub = "data_rec"

        return BASE_PREP_DIR / sub

    # everything else uses DATA_DIR_MAP
    return Path(DATA_DIR_MAP[tag])


def screen_ls():
    try:
        out = subprocess.check_output(
            ["screen", "-ls"],
            text=True,
            stderr=subprocess.STDOUT,
        )
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


def start_screen_session(name: str,
                         cmd: str,
                         log_path: Path = None,
                         workdir: Path = None):
    """
    Launch a detached screen session running `cmd`.
    The session will EXIT when `cmd` finishes because we use `exec`.
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
        choices=[
            "Data", "Rec", "Gen",
            "DATA", "REC", "GEN",
            "Dis_data", "Dis_rec", "DIS_GEN",
        ],
        default="Data",
        help="Which dataset to process: Data, Rec, or Gen (case-insensitive).",
    )

    parser.add_argument(
        "--is-true-gen-event",
        type=int,
        choices=[0, 1],
        default=1,
        help="Pass 0 or 1 to the macro as the is_true_gen_event flag.",
    )

    parser.add_argument(
        "--mx-cut-mode",
        type=int,
        choices=[0, 1, 2],
        default=2,
        help="Mx cut mode: 0 = no cut, 1 = Mx > 1.0, 2 = Mx > 1.5.",
    )

    parser.add_argument(
        "--custom-output",
        type=int,
        choices=[0, 1],
        default=0,
        help="Use data_test input directory and save custom output columns in the Data snapshot.",
    )

    parser.add_argument(
        "--isGoodElec",
        type=int,
        choices=[0, 1],
        default=0,
        help="Requires detected electron for SIDIS acceptance.",
    )

    parser.add_argument(
        "--isMCParentCut",
        type=int,
        choices=[0, 1],
        default=0,
        help="Gammas have to originate from pi0, MC only.",
    )

    parser.add_argument(
        "--noPhiBinning",
        type=int,
        choices=[0, 1],
        default=0,
        help="Integrate over phi.",
    )

    parser.add_argument(
        "--gregs_ai_on",
        type=int,
        choices=[0, 1],
        default=0,
        help="If AI ID was used, affects cuts and output path.",
    )

    parser.add_argument(
        "--gg_mom_cut",
        type=int,
        choices=[0, 1],
        default=0,
        help="0 is g_mom_cut > 0.5 GeV, 1 is > 0.35 GeV. Default is 0.5 GeV cut.",
    )

    parser.add_argument(
        "--no_cuts_gen",
        type=int,
        choices=[0, 1],
        default=0,
        help="For acceptance studies. Does not apply DIS cuts to SIDIS but saves REC electron information with the cuts columns.",
    )

    parser.add_argument(
        "--logs-dir",
        default="logs_add_binning",
        help="Directory to store screen logs.",
    )

    args = parser.parse_args()

    data_tag = args.data_type.upper()

    if data_tag not in ("DATA", "REC") and data_tag not in DATA_DIR_MAP:
        raise SystemExit(f"[ERROR] Unknown data type '{args.data_type}'.")

    data_dir_path = pick_input_dir(
        data_tag=data_tag,
        mx_cut_mode=int(args.mx_cut_mode),
        is_true_gen_event=int(args.is_true_gen_event),
        is_mc_parent_cuts=int(args.isMCParentCut),
        no_phi_binning=int(args.noPhiBinning),
        gregs_ai_on=int(args.gregs_ai_on),
        gg_mom_cut=int(args.gg_mom_cut),
        custom_output=int(args.custom_output),
    )

    # Validate input dir
    mypath = Path(data_dir_path)
    print(data_dir_path)

    if not mypath.is_dir():
        raise SystemExit(f"[ERROR] Input directory not found: {mypath}")

    # Gather input files
    rec_files = [p.name for p in mypath.iterdir() if p.is_file()]

    if not rec_files:
        raise SystemExit(f"[ERROR] No input files found in: {mypath}")

    logs_dir = Path(args.logs_dir)

    is_true_gen_event = int(args.is_true_gen_event)
    mx_cut_mode = int(args.mx_cut_mode)
    custom_output = int(args.custom_output)
    isGoodElec = int(args.isGoodElec)
    isMCParentCut = int(args.isMCParentCut)
    noPhiBinning = int(args.noPhiBinning)
    gregs_ai_on = int(args.gregs_ai_on)
    gg_mom_cut = int(args.gg_mom_cut)
    no_cuts_gen = int(args.no_cuts_gen)

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
        name = f"{SCREEN_NAME}_{i_file}"
        log = logs_dir / f"{name}.log"

        # Build the macro call.
        # Signature:
        # add_binning.cxx(
        #   "<file>",
        #   "<DATA|REC|GEN>",
        #   is_true_gen_event,
        #   mx_cut_mode,
        #   custom_output,
        #   isGoodElec,
        #   isMCParentCut,
        #   noPhiBinning,
        #   gregs_ai_on,
        #   gg_mom_cut,
        #   no_cuts_gen
        # )
        macro_call = (
            f'source/add_binning.cxx('
            f'"{file_path}", '
            f'"{data_tag}", '
            f'{is_true_gen_event}, '
            f'{mx_cut_mode}, '
            f'{custom_output}, '
            f'{isGoodElec}, '
            f'{isMCParentCut}, '
            f'{noPhiBinning}, '
            f'{gregs_ai_on}, '
            f'{gg_mom_cut}, '
            f'{no_cuts_gen})'
        )

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

        print(
            f"Launched {name} "
            f"[data_type={data_tag}] "
            f"(is_true_gen_event={is_true_gen_event}) "
            f"(custom_output={custom_output}) "
            f"(log: {log})"
        )