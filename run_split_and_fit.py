#!/usr/bin/env python3
import argparse
import os
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path

MAX_SCREENS = 25
poll_sec = 1
SCREEN_NAME = "split_fit"

BASE_PATH = Path("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/")

def parse_args():
    ap = argparse.ArgumentParser(description="Launch split_and_fit ROOT jobs in screen sessions.")
    ap.add_argument(
        "--subdir",
        choices=["unfolding_rec_true", "unfolding_rec_fake"],
        default="unfolding_rec_true",
        help="Folder to append to base path for input files."
    )
    return ap.parse_args()

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
    args = parse_args()

    LD_EXPORT = (
      'export LD_LIBRARY_PATH='
      '"/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.6b/4.3.0/lib64:'
      '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.7/lib:'
      '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.36.04/lib:'
      '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/python/3.13.7/lib:'
      '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.3.0/lib:'
      '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib64:'
      '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib:${LD_LIBRARY_PATH-}"'
    )

    # Append requested subfolder to base path
    mypath = BASE_PATH / args.subdir
    if not mypath.exists():
        print(f"WARNING: input directory does not exist: {mypath}", flush=True)

    # match a minus sign directly followed by digits, anywhere in the name
    has_negative = re.compile(r'-(?=\d)')

    rec_files = [
        p.name
        for p in mypath.iterdir() if mypath.exists()
        if p.is_file()
        and "fitted" not in p.name.lower()      # exclude anything with "fitted"
        and not has_negative.search(p.stem)     # skip names containing a negative number
    ]

    logs_dir = Path("logs_split_and_fit")

    for i_file, file in enumerate(rec_files):
        # wait until there is available screen slot:
        while len(screen_ls()) >= MAX_SCREENS:
            time.sleep(poll_sec)

        file_path = str(mypath / file)
        name = f"{SCREEN_NAME}_{i_file}"        # screen session name must be unique
        log = logs_dir / f"{name}.log"          # capture stdout/stderr

        # Escape any embedded quotes in the file path for the ROOT macro call
        file_for_root = file_path.replace('"', '\\"')

        # Build the macro call and ROOT command
        macro_call = f'source/split_and_fit.cxx("{file_for_root}")'
        root_call = f"root -l -b -q {shlex.quote(macro_call)}"

        cmd = f"""
        echo "Start: $(date)"
        {LD_EXPORT}
        {root_call}
        echo "Done: $(date)"
        """

        start_screen_session(name, cmd, log_path=log)
        print(f"Launched {name}: {cmd}  (log: {log})")
