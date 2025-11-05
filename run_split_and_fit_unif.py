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
    '"/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.6b/4.3.0/lib64:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.7/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.36.04/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/python/3.13.7/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.3.0/lib:'
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
    ap.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parent,  # folder that contains ./source/
        help="Working directory for the screen job; should contain 'source/'.",
    )
    ap.add_argument(
        "--aclic-dir",
        type=Path,
        default=Path("aclic_build"),
        help="Where to put ACLiC outputs (per job subfolder will be created here).",
    )
    ap.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Use '++' on .L to force rebuild of the macro.",
    )
    return ap.parse_args()

def screen_ls():
    try:
        out = subprocess.check_output(["screen", "-ls"], text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        out = e.output

    sessions = {}
    for line in out.splitlines():
        # lines like: "    12345.name   (Detached)"
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

    project_root = args.project_root.resolve()
    macro_path = (project_root / "source" / "split_and_fit_unified.cxx").resolve()

    for i_file, file in enumerate(sorted(rec_files)):
        # throttle concurrent sessions
        while len(screen_ls()) >= args.max_screens:
            time.sleep(POLL_SEC)

        file_path = str((mypath / file).resolve())
        name = f"{SCREEN_NAME}_{i_file:04d}"     # unique screen session
        log = LOGS_DIR / f"{name}.log"           # capture stdout/stderr

        # Escape any embedded quotes in the file path for the ROOT call
        file_for_root = file_path.replace('"', '\\"')

        # Unique ACLiC build dir per job to avoid collisions
        build_dir = (project_root / args.aclic_dir / name).resolve()

        # ROOT command lines
        build_line = f'gSystem->SetBuildDir("{build_dir}", kTRUE)'
        plus = "++" if args.force_rebuild else "+"
        load_line = f'.L {macro_path}{plus}'
        run_line  = f'split_and_fit_switch("{file_for_root}", "{logic}", {int(args.png_every)})'

        # Ensure build dir exists before ROOT starts
        mkdir_cmd = f'mkdir -p {shlex.quote(str(build_dir))}'

        root_call = (
            f"root -l -b -q "
            f"-e {shlex.quote(build_line)} "
            f"-e {shlex.quote(load_line)} "
            f"-e {shlex.quote(run_line)}"
        )

        cmd = f"""
        echo "Start: $(date)"
        {LD_EXPORT}
        {mkdir_cmd}
        {root_call}
        echo "Done: $(date)"
        """

        # Run with project_root as working directory so relative paths (e.g. 'source/') are stable
        start_screen_session(name, cmd, log_path=log, workdir=project_root)
        print(f"Launched {name}: (cwd={project_root})  (log: {log})")
