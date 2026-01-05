#!/usr/bin/env python3
# Launch clas12root gen-dis jobs in detached screen sessions with module setup and throttling.
# Scans one or more base directories recursively for *.hipo files.

import sys, re, time, shutil, subprocess
from pathlib import Path
from datetime import datetime

# --- settings (edit as needed) ---
BASE_DIRS = [
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/",
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/Q2_1.5GeV/"
]

FILE_GLOBS = ("*.hipo",)

# DIS-specific output
rootpath   = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec_dis/"
rootprefix = "rec_dis_f2018_"

# ROOT macro for DIS generator output
macro = "source/gen_root_files_rec_dis.C"

MAX_PARALLEL = 25

# --- basic checks ---
if shutil.which("screen") is None:
    sys.exit("ERROR: 'screen' not found in PATH.")

clean_bases: list[Path] = []
for d in BASE_DIRS:
    p = Path(d).resolve()
    if not p.is_dir():
        sys.exit(f"ERROR: Base directory not found: {p}")
    clean_bases.append(p)

# NEW: ensure children (like .../fa18_inb/Q2_1.5GeV) are scanned before parents (.../fa18_inb)
clean_bases.sort(key=lambda p: len(p.parts), reverse=True)

if not rootpath.endswith("/"):
    rootpath += "/"

# --- scan files recursively in all trees; dedupe by real path ---
def iter_files(base: Path):
    for pat in FILE_GLOBS:
        for f in base.rglob(pat):
            if f.is_file():
                rp = f.resolve()
                rel = rp.relative_to(base)
                yield base, rp, rel

seen = set()
files: list[tuple[Path, Path, Path]] = []
for base in clean_bases:
    for base_dir, abs_path, rel_path in iter_files(base):
        if abs_path in seen:
            continue
        seen.add(abs_path)
        files.append((base_dir, abs_path, rel_path))

if not files:
    sys.exit("Nothing to do (no matching files found under BASE_DIRS).")

print(f"Found {len(files)} input .hipo files under:")
for b in clean_bases:
    print("  ", b)

# --- setup logging & job prefix ---
JOB_PREFIX = f"gen_dis_{datetime.now():%Y%m%d_%H%M%S}"
LOG_DIR    = Path("logs"); LOG_DIR.mkdir(exist_ok=True)
MISSING    = Path("gen_dis_file.log");        MISSING.write_text("")
LAUNCHED   = Path("launched_gen_dis_jobs.log"); LAUNCHED.write_text("")

def running_jobs() -> int:
    """Count running screen sessions for this batch."""
    out = subprocess.run(
        ["screen", "-ls"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True
    )
    return out.stdout.count(JOB_PREFIX)

# env + modules run inside each screen (same style as your newer script)
module_cmd = (
    '[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh; '
    'module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles; '
    'module load clas12root; module load root'
)

# LD_LIBRARY_PATH inside each screen (updated to new versions)
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

# --- main loop over discovered files ---
for base_dir, abs_path, rel_path in files:
    hipopath_arg = (str(base_dir) if str(base_dir).endswith("/") else str(base_dir) + "/")
    f_rel        = str(rel_path)

    if not abs_path.is_file():
        print(f"Missing: {abs_path}")
        MISSING.write_text(MISSING.read_text() + str(abs_path) + "\n")
        continue

    print(f"Found: {abs_path}")

    # ROOT command: gen_root_files_gen_dis.C(hipopath, relative file, rootpath+rootprefix)
    cmd  = (
        "clas12root -q "
        f"'{macro}(\"{hipopath_arg}\", \"{f_rel}\", \"{rootpath}{rootprefix}\")'"
    )

    # Screen session name + logfile
    safe = re.sub(r"[^A-Za-z0-9_-]", "_", rel_path.name)
    sess = f"{JOB_PREFIX}_{safe}"
    log  = LOG_DIR / f"{sess}.log"

    # run modules -> export LD_LIBRARY_PATH -> run ROOT macro
    full = f'{module_cmd}; {LD_EXPORT}; {cmd} > "{log}" 2>&1'
    subprocess.run(["screen", "-dmS", sess, "bash", "-lc", full], check=True)

    with LAUNCHED.open("a") as lf:
        lf.write(f"{sess}\tbase={base_dir}\tfile={f_rel}\tlog={log}\n")

    # throttle to MAX_PARALLEL
    while running_jobs() >= MAX_PARALLEL:
        time.sleep(2)

print()
print(f"Launched screen sessions with prefix: {JOB_PREFIX}")
print(f"Missing files -> {MISSING}")
print(f"Launched jobs log -> {LAUNCHED}")
print(f"Logs -> {LOG_DIR}")
print()
print("Handy commands:")
print(f"  screen -ls | grep {JOB_PREFIX}")
print("  screen -r <session_name>")
print("  screen -S <session_name> -X quit")
print(f"  tail -f {LOG_DIR}/{JOB_PREFIX}_<sanitized-filename>.log")
