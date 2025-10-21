#!/usr/bin/env python3
import argparse
import re
import shlex
import subprocess
import time
from pathlib import Path

# ---------------- config ----------------
MAX_SCREENS = 25
POLL_SEC = 1
SCREEN_NAME = "split_fit"

BASE_PATH = Path("/w/hallb-scshelf2102/clas12/valerii/multiPi0/pass2_v3/")
LOGS_DIR = Path("logs_split_and_fit")

# ROOT env (same as your originals)
LD_EXPORT = (
    'export LD_LIBRARY_PATH='
    '"/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.30.04/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.5/4.2.0/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/iguana/0.8.0/4.2.0/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib64:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib:${LD_LIBRARY_PATH-}"'
)

# ---------------- helpers ----------------
def parse_args():
    ap = argparse.ArgumentParser(
        description="Launch unified split_and_fit ROOT jobs in screen sessions."
    )
    ap.add_argument(
        "--subdir",
        choices=["unfolding_rec_data", "unfolding_rec_true", "unfolding_rec_fake"],
        default="unfolding_rec_data",
        help="Folder (under BASE_PATH) containing input .root files.",
    )
    ap.add_argument(
        "--logic",
        choices=["auto", "data", "sim"],
        default="auto",
        help="Which unified logic to use. Default 'auto' infers from --subdir.",
    )
    ap.add_argument(
        "--png-every",
        type=int,
        default=-1,
        help="Pass-through to C++: save every Nth PNG (<=0 keeps C++ defaults: 50 data, 2 sim).",
    )
    ap.add_argument(
        "--max-screens",
        type=int,
        default=MAX_SCREENS,
        help="Maximum concurrent screen sessions.",
    )
    return ap.parse_args()

def screen_ls():
    try:
        out = subprocess.check_output(["screen", "-ls"], text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        out = e.output

    sessions = {}
    for line in out.splitlines():
        # lines look like: "    12345.name   (Detached)"
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

# ---------------- main ----------------
if __name__ == "__main__":
    args = parse_args()

    # infer logic from subdir if requested
    logic = args.logic
    if logic == "auto":
        logic = "data" if args.subdir == "unfolding_rec_data" else "sim"

    # build input directory
    mypath = BASE_PATH / args.subdir
    if not mypath.exists():
        print(f"WARNING: input directory does not exist: {mypath}", flush=True)

    # filter: exclude "fitted" and any name containing a negative number
    has_negative = re.compile(r'-(?=\d)')

    rec_files = []
    if mypath.exists():
        for p in mypath.iterdir():
            if p.is_file() and ("fitted" not in p.name.lower()) and (not has_negative.search(p.stem)):
                rec_files.append(p.name)

    for i_file, file in enumerate(sorted(rec_files)):
        # throttle concurrent sessions
        while len(screen_ls()) >= args.max_screens:
            time.sleep(POLL_SEC)

        file_path = str(mypath / file)
        name = f"{SCREEN_NAME}_{i_file:04d}"     # unique screen session
        log = LOGS_DIR / f"{name}.log"           # capture stdout/stderr

        # Escape any embedded quotes in the file path for the ROOT call
        file_for_root = file_path.replace('"', '\\"')

        # We load the unified macro once then call the switch entry point.
        #   split_and_fit_switch(path, "data"|"sim", png_every)
        macro_cmd = (
            f'source/split_and_fit_unified.cxx; '
            f'split_and_fit_switch("{file_for_root}", "{logic}", {int(args.png_every)})'
        )
        root_call = f"root -l -b -q {shlex.quote(macro_cmd)}"

        cmd = f"""
        echo "Start: $(date)"
        {LD_EXPORT}
        {root_call}
        echo "Done: $(date)"
        """

        start_screen_session(name, cmd, log_path=log)
        print(f"Launched {name}: {cmd}  (log: {log})")
