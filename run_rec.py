#!/usr/bin/env python3
import os, sys, shutil, subprocess
from pathlib import Path
from datetime import datetime

# --- user settings ---
BASE_DIRS = [
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/",
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/Q2_1.5GeV/"
]
# file patterns to include (adjust as needed)
FILE_GLOBS = ("*.hipo")

rootpath   = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/"
rootprefix = "rec_f2018_"
MAX_PARALLEL = 35

# --- checks ---
if shutil.which("screen") is None:
    sys.exit("ERROR: 'screen' not found in PATH.")

# validate base dirs
clean_base_dirs = []
for d in BASE_DIRS:
    p = Path(d).resolve()
    if not p.is_dir():
        sys.exit(f"ERROR: Base directory not found: {p}")
    clean_base_dirs.append(p)

# --- scan files recursively in both trees ---
def iter_files(base: Path):
    seen = set()
    for pat in FILE_GLOBS:
        for f in base.rglob(pat):
            if f.is_file():
                # dedupe by real path
                rp = f.resolve()
                if rp in seen:
                    continue
                seen.add(rp)
                rel = rp.relative_to(base)
                yield {
                    "base": base,                    # base dir
                    "rel": str(rel),                 # relative path under base
                    "name": rp.name,                 # filename
                    "stem": rp.stem,                 # filename without extension
                    "abs": rp                        # absolute path
                }

files = []
for base in clean_base_dirs:
    files.extend(iter_files(base))

if not files:
    sys.exit("Nothing to do (no matching files found).")

# --- setup ---
JOB_PREFIX   = f"rec_{datetime.now():%Y%m%d_%H%M%S}"
LOG_DIR      = Path("logs").resolve(); LOG_DIR.mkdir(exist_ok=True)
BATCH_DIR    = Path(f"batches_{JOB_PREFIX}").resolve(); BATCH_DIR.mkdir(exist_ok=True)
MISSING_LOG  = Path("rec_file.log").resolve();  MISSING_LOG.write_text("")
LAUNCHED_LOG = Path("launched_rec_jobs.log").resolve(); LAUNCHED_LOG.write_text("")

# bucketize work
BATCH_COUNT = min(MAX_PARALLEL, len(files))
buckets = [[] for _ in range(BATCH_COUNT)]
for i, meta in enumerate(files):
    # TSV: base<TAB>rel<TAB>name<TAB>stem
    line = f"{meta['base']}\t{meta['rel']}\t{meta['name']}\t{meta['stem']}"
    buckets[i % BATCH_COUNT].append(line)

# LD_LIBRARY_PATH inside each screen (keep if you still want it alongside modules)
LD_EXPORT = (
    'export LD_LIBRARY_PATH="/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.30.04/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.5/4.2.0/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/iguana/0.8.0/4.2.0/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib64:'
    '/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib:${LD_LIBRARY_PATH-}"'
)

for b, bucket in enumerate(buckets):
    if not bucket: continue
    lst = BATCH_DIR / f"batch_{b}.tsv"
    lst.write_text("\n".join(bucket) + "\n")

    runner = BATCH_DIR / f"runner_{JOB_PREFIX}_b{b}.sh"
    runner.write_text(f"""#!/usr/bin/env bash
set -euo pipefail

# --- environment modules first ---
[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12root
module load root

JOB_PREFIX="{JOB_PREFIX}"
ROOTPATH="{rootpath}"
ROOTPREFIX="{rootprefix}"
LOG_DIR="{LOG_DIR}"
MISSING_LOG="{MISSING_LOG}"
BATCH_LIST="{lst}"

{LD_EXPORT}

while IFS=$'\\t' read -r HIPOPATH file name stem || [[ -n "$HIPOPATH" ]]; do
  [[ -z "${{HIPOPATH:-}}" || -z "${{file:-}}" ]] && continue

  
  fpath="${{HIPOPATH%/}}/${{file}}"
  safe_rel=$(echo "$file" | tr -c 'A-Za-z0-9_./-' '_')
  safe=$(basename "$safe_rel" | tr -c 'A-Za-z0-9_-' '_')
  LOG_FILE="${{LOG_DIR}}/${{JOB_PREFIX}}_${{safe}}.log"

  if [[ ! -f "$fpath" ]]; then
    echo "$fpath" >> "$MISSING_LOG"
    continue
  fi

  # If your macro only takes 3 args, use this instead:
  cmd="clas12root -q 'source/gen_root_files_rec.C(\\"${{HIPOPATH}}\\", \\"${{file}}\\", \\"${{ROOTPATH}}${{ROOTPREFIX}}\\")'"

  bash -lc "$cmd" > "$LOG_FILE" 2>&1
done < "$BATCH_LIST"
""")
    runner.chmod(0o755)

    sess = f"{JOB_PREFIX}_batch{b}"
    batch_log = LOG_DIR / f"{sess}.log"
    subprocess.run(["screen", "-dmS", sess, "bash", "-lc", f'"{runner}" > "{batch_log}" 2>&1'], check=True)
    with LAUNCHED_LOG.open("a") as lf:
        lf.write(f"{sess}\tbatch={b}\tfiles={len(bucket)}\trunner={runner}\tlog={batch_log}\n")

print(f"Launched {BATCH_COUNT} screen session(s) with prefix: {JOB_PREFIX}")
print(f"Missing files: {MISSING_LOG}")
print(f"Launched jobs: {LAUNCHED_LOG}")
print(f"Per-file logs in: {LOG_DIR}")
