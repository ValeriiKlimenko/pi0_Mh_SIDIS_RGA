#!/usr/bin/env python3
# Launch clas12root gen jobs in detached screen sessions with module setup and throttling.
# Scans only the selected base directory (non-recursive) for *.hipo files.

import sys, re, time, shutil, subprocess, argparse
from pathlib import Path
from datetime import datetime

# --- settings ---
BASE_DIRS = [
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/",
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/Q2_1.5GeV/",
]

FILE_GLOBS = ("*.hipo",)

rootpath   = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/"
rootprefix = "gen_f2018_"
MAX_PARALLEL = 10

# --- CLI ---
ap = argparse.ArgumentParser(
    description="Launch clas12root gen jobs from one selected BASE_DIR entry, ignoring subfolders."
)
ap.add_argument(
    "--base-index",
    type=int,
    default=0,
    help="Which entry of BASE_DIRS to use (default: 0)."
)
args = ap.parse_args()

if not (0 <= args.base_index < len(BASE_DIRS)):
    sys.exit(
        f"ERROR: --base-index must be between 0 and {len(BASE_DIRS)-1}. "
        f"Got {args.base_index}."
    )

selected_base = Path(BASE_DIRS[args.base_index]).resolve()

# --- checks ---
if shutil.which("screen") is None:
    sys.exit("ERROR: 'screen' not found in PATH.")

if not selected_base.is_dir():
    sys.exit(f"ERROR: Base directory not found: {selected_base}")

# --- scan files only in selected base dir, not in subfolders ---
def iter_files(base: Path):
    for pat in FILE_GLOBS:
        for f in base.glob(pat):   # non-recursive
            if f.is_file() and f.parent.resolve() == base:
                rp = f.resolve()
                rel = Path(f.name)  # just the filename, no subdir component
                yield base, rp, rel

files: list[tuple[Path, Path, Path]] = list(iter_files(selected_base))

if not files:
    sys.exit(f"Nothing to do (no matching files found in top level of {selected_base}).")

# --- setup ---
JOB_PREFIX = f"gen_{datetime.now():%Y%m%d_%H%M%S}"
LOG_DIR    = Path("logs")
LOG_DIR.mkdir(exist_ok=True)

MISSING  = Path("gen_file.log")
MISSING.write_text("")

LAUNCHED = Path("launched_rec_jobs.log")
LAUNCHED.write_text("")

def running_jobs() -> int:
    out = subprocess.run(
        ["screen", "-ls"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True
    )
    return out.stdout.count(JOB_PREFIX)

# env + modules run inside each screen
module_cmd = (
    '[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh; '
    'module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles; '
    'module load clas12root; module load root'
)

# LD_LIBRARY_PATH inside each screen
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

for base_dir, abs_path, rel_path in files:
    hipopath_arg = str(base_dir)
    if not hipopath_arg.endswith("/"):
        hipopath_arg += "/"

    f_rel = str(rel_path)   # filename only

    if not abs_path.is_file():
        with MISSING.open("a") as mf:
            mf.write(str(abs_path) + "\n")
        continue

    safe = re.sub(r"[^A-Za-z0-9_-]", "_", rel_path.name)
    sess = f"{JOB_PREFIX}_{safe}"
    log  = LOG_DIR / f"{sess}.log"

    # pass: (hipopath, relative file under that base, output prefix)
    cmd = (
        "clas12root -q "
        f"'source/gen_root_files_gen.C(\"{hipopath_arg}\", \"{f_rel}\", \"{rootpath}{rootprefix}\")'"
    )

    full = f'{module_cmd}; {LD_EXPORT}; {cmd} > "{log}" 2>&1'
    subprocess.run(["screen", "-dmS", sess, "bash", "-lc", full], check=True)

    with LAUNCHED.open("a") as lf:
        lf.write(f"{sess}\tbase={base_dir}\tfile={f_rel}\tlog={log}\n")

    while running_jobs() >= MAX_PARALLEL:
        time.sleep(2)

print(f"Using BASE_DIRS[{args.base_index}] = {selected_base}")
print(f"Launched jobs with prefix: {JOB_PREFIX}")
print(f"Missing files -> {MISSING}")
print(f"Launched jobs log -> {LAUNCHED}")
print(f"Logs -> {LOG_DIR}")