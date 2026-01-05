#!/usr/bin/env python3
import os, sys, shutil, subprocess, socket, getpass
from pathlib import Path
from datetime import datetime

# --- user settings ---
BASE_DIRS = [
    #"/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/",
    "/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/Q2_1.5GeV/"
]
# file patterns to include (adjust as needed)
FILE_GLOBS = ("*.hipo",)

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
                rp = f.resolve()
                if rp in seen:
                    continue
                seen.add(rp)
                rel = rp.relative_to(base)
                yield {
                    "base": base,
                    "rel": str(rel),
                    "name": rp.name,
                    "stem": rp.stem,
                    "abs": rp
                }

files = []
per_base_counts = {}
for base in clean_base_dirs:
    items = list(iter_files(base))
    files.extend(items)
    per_base_counts[str(base)] = len(items)

if not files:
    sys.exit("Nothing to do (no matching files found).")

# --- setup ---
JOB_PREFIX   = f"rec_{datetime.now():%Y%m%d_%H%M%S}"
LOG_DIR      = Path("logs").resolve(); LOG_DIR.mkdir(exist_ok=True)
BATCH_DIR    = Path(f"batches_{JOB_PREFIX}").resolve(); BATCH_DIR.mkdir(exist_ok=True)
MISSING_LOG  = Path("rec_file.log").resolve();  MISSING_LOG.write_text("")
FAILED_LOG   = Path("failed_rec_jobs.log").resolve(); FAILED_LOG.write_text("")
LAUNCHED_LOG = Path("launched_rec_jobs.log").resolve(); LAUNCHED_LOG.write_text("")

# --- debug banner ---
print("="*78)
print(f"[{datetime.now():%F %T}] Launcher host={socket.gethostname()} user={getpass.getuser()}")
print(f"Working dir: {Path.cwd()}")
print(f"JOB_PREFIX:  {JOB_PREFIX}")
print("Base directories & counts:")
for b, n in per_base_counts.items():
    print(f"  - {b} : {n} files")
print(f"Total files found: {len(files)}")
print(f"LOG_DIR:   {LOG_DIR}")
print(f"BATCH_DIR: {BATCH_DIR}")
print("="*78)

# bucketize work
BATCH_COUNT = min(MAX_PARALLEL, len(files))
buckets = [[] for _ in range(BATCH_COUNT)]
for i, meta in enumerate(files):
    # TSV: base<TAB>rel<TAB>name<TAB>stem
    line = f"{meta['base']}\t{meta['rel']}\t{meta['name']}\t{meta['stem']}"
    buckets[i % BATCH_COUNT].append(line)

# LD_LIBRARY_PATH inside each screen (keep if you still want it alongside modules)
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

for b, bucket in enumerate(buckets):
    if not bucket: 
        print(f"[launcher] bucket {b} is empty, skipping.")
        continue

    lst = BATCH_DIR / f"batch_{b}.tsv"
    lst.write_text("\n".join(bucket) + "\n")
    print(f"[launcher] wrote {lst} with {len(bucket)} lines")

    runner = BATCH_DIR / f"runner_{JOB_PREFIX}_b{b}.sh"
    runner.write_text(f"""#!/usr/bin/env bash
set -euo pipefail

START_TS=$(date +%s)
shopt -s lastpipe

# --- environment modules first ---
[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh || true

JOB_PREFIX="{JOB_PREFIX}"
ROOTPATH="{rootpath}"
ROOTPREFIX="{rootprefix}"
LOG_DIR="{LOG_DIR}"
MISSING_LOG="{MISSING_LOG}"
FAILED_LOG="{FAILED_LOG}"
BATCH_LIST="{lst}"

{LD_EXPORT}

# --- runner banner (goes to batch log) ---
echo "======================================================================="
echo "[runner] START $(date '+%F %T')"
echo "[runner] host=$(hostname) user=$(whoami) pid=$$"
echo "[runner] pwd=$PWD"
echo "[runner] PATH=$PATH"
echo "[runner] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "[runner] BATCH_LIST=$BATCH_LIST"
echo "[runner] LOG_DIR=$LOG_DIR"
echo "[runner] ROOTPATH=$ROOTPATH ROOTPREFIX=$ROOTPREFIX"
echo "======================================================================="

# sanity: clas12root available?
if ! command -v clas12root >/dev/null 2>&1; then
  echo "[runner][FATAL] 'clas12root' not found in PATH" >&2
  exit 127
fi
echo "[runner] clas12root=$(command -v clas12root)"
# try to report a version (non-fatal if it errors)
clas12root -q -b -l -x 'gSystem->ProcessEvents(); gApplication->Terminate(0);' >/dev/null 2>&1 || true

TOTAL=$(wc -l < "$BATCH_LIST" || echo 0)
i=0
ok=0
fail=0
missing=0

# loop
while IFS=$'\\t' read -r HIPOPATH file name stem || [[ -n "${{HIPOPATH:-}}" ]]; do
  [[ -z "${{HIPOPATH:-}}" || -z "${{file:-}}" ]] && continue
  i=$((i+1))

  fpath="${{HIPOPATH%/}}/${{file}}"
  safe_rel=$(echo "$file" | tr -c 'A-Za-z0-9_./-' '_')
  safe=$(basename "$safe_rel" | tr -c 'A-Za-z0-9_-' '_')
  LOG_FILE="${{LOG_DIR}}/${{JOB_PREFIX}}_${{safe}}.log"

  echo "[{b}][$i/$TOTAL] $(date '+%F %T') file=$file abs=$fpath"

  if [[ ! -f "$fpath" ]]; then
    echo "$fpath" >> "$MISSING_LOG"
    echo "[{b}][$i/$TOTAL] MISSING: $fpath (logged to $MISSING_LOG)"
    missing=$((missing+1))
    continue
  fi

  # If your macro only takes 3 args, use this instead:
  cmd="clas12root -q 'source/gen_root_files_rec.C(\\"${{HIPOPATH}}\\", \\"${{file}}\\", \\"${{ROOTPATH}}${{ROOTPREFIX}}\\")'"

  echo "[{b}][$i/$TOTAL] START $(date '+%F %T') → $cmd" | tee -a "$LOG_FILE"
  SECONDS=0
  set +e
  bash -lc "$cmd" >> "$LOG_FILE" 2>&1
  rc=$?
  set -e
  dur=$SECONDS

  if [[ $rc -eq 0 ]]; then
    echo "[{b}][$i/$TOTAL] DONE rc=$rc time=${{dur}}s" | tee -a "$LOG_FILE"
    ok=$((ok+1))
  else
    echo "[{b}][$i/$TOTAL] FAIL rc=$rc time=${{dur}}s (see $LOG_FILE)" | tee -a "$LOG_FILE"
    echo -e "$rc\\t$fpath\\t$LOG_FILE" >> "$FAILED_LOG"
    fail=$((fail+1))
  fi
done < "$BATCH_LIST"

END_TS=$(date +%s)
echo "======================================================================="
echo "[runner] END   $(date '+%F %T')  ok=$ok fail=$fail missing=$missing total=$TOTAL elapsed=$((END_TS-START_TS))s"
echo "[runner] missing log: $MISSING_LOG"
echo "[runner] failed  log: $FAILED_LOG"
echo "======================================================================="
""")
    runner.chmod(0o755)

    sess = f"{JOB_PREFIX}_batch{b}"
    batch_log = LOG_DIR / f"{sess}.log"
    # note: redirect happens inside the shell passed to screen
    subprocess.run(
        ["screen", "-dmS", sess, "bash", "-lc", f'"{runner}" > "{batch_log}" 2>&1'],
        check=True
    )
    with LAUNCHED_LOG.open("a") as lf:
        lf.write(f"{sess}\tbatch={b}\tfiles={len(bucket)}\trunner={runner}\tlog={batch_log}\n")
    print(f"[launcher] launched session={sess} files={len(bucket)} runner={runner}")
    print(f"           batch log → {batch_log}")

print("="*78)
print(f"Launched {sum(1 for buk in buckets if buk)} screen session(s) with prefix: {JOB_PREFIX}")
print(f"Missing files log : {MISSING_LOG}")
print(f"Failed jobs log   : {FAILED_LOG}")
print(f"Launched jobs list: {LAUNCHED_LOG}")
print(f"Per-file logs in  : {LOG_DIR}")
print("Tips:")
print("  • screen -ls")
print(f"  • screen -r {JOB_PREFIX}_batch0   # attach to first batch")
print("  • tail -f logs/<session>.log       # live batch runner output")
print("  • grep -v '^#' failed_rec_jobs.log | cut -f3 | xargs -I{{}} tail -n 3 {{}}")
print("="*78)
