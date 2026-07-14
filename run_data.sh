#!/usr/bin/env bash
# Run gen_root_files_data.C jobs in parallel inside detached GNU screen sessions
# and write each session's output to logs/<session>.log

set -euo pipefail

# -----------------------------
# CLI option: --is_greg_ai 0/1
# -----------------------------
IS_GREG_AI=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --is_greg_ai)
      [[ $# -ge 2 ]] || { echo "ERROR: --is_greg_ai needs 0 or 1" >&2; exit 2; }
      IS_GREG_AI="$2"
      shift 2
      ;;
    --is_greg_ai=*)
      IS_GREG_AI="${1#*=}"
      shift 1
      ;;
    -h|--help)
      cat <<EOF
Usage: $0 [--is_greg_ai 0|1]

  --is_greg_ai 0   (default) output to .../data_data/
  --is_greg_ai 1             output to .../data_data_greg_ai/
EOF
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ "$IS_GREG_AI" != "0" && "$IS_GREG_AI" != "1" ]]; then
  echo "ERROR: --is_greg_ai must be 0 or 1 (got '$IS_GREG_AI')" >&2
  exit 2
fi

# --- user settings (kept from your script) ---
input="path_to_files/pass2_nSidis.dat"
hipopath="/lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/"

# rootpath depends on is_greg_ai
if [[ "$IS_GREG_AI" == "1" ]]; then
  rootpath="/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data_greg_ai/"
else
  rootpath="/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/"
fi

mkdir -p "$rootpath"

rootprefix="data_f2018_"

# How many screen sessions to run at once
MAX_PARALLEL=30

# --- safety checks & prep ---
command -v screen >/dev/null 2>&1 || { echo "ERROR: 'screen' is not installed/in PATH." >&2; exit 1; }
[[ -f "$input" ]] || { echo "ERROR: Input file list not found: $input" >&2; exit 1; }

JOB_PREFIX="nSidis_$(date +%Y%m%d_%H%M%S)"
MISSING_LOG="data_file.log"        # matches your original name
LAUNCHED_LOG="launched_jobs.log"
LOG_DIR="logs"

mkdir -p "$LOG_DIR"
: > "$MISSING_LOG"
: > "$LAUNCHED_LOG"

echo "[launcher] is_greg_ai=$IS_GREG_AI"
echo "[launcher] rootpath=$rootpath"

# Load file list
mapfile -t files < "$input"
echo "Read ${#files[@]} entries from $input"

# Helper: count running screen jobs for this batch
running_jobs() {
  screen -ls 2>/dev/null | grep -c "$JOB_PREFIX" || true
}

# LD_LIBRARY_PATH export to run inside each screen (note the escaped quotes and \$)
LD_EXPORT="export LD_LIBRARY_PATH=\"/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.9.0/4.3.0/lib64:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.7/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.36.04/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/python/3.13.7/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.3.0/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib64:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib:\${LD_LIBRARY_PATH-}\""

# --- main loop ---
for file in "${files[@]}"; do
  # skip empty lines / comments
  [[ -z "$file" || "$file" =~ ^[[:space:]]*# ]] && continue

  fpath="${hipopath}${file}"
  if [[ ! -f "$fpath" ]]; then
    echo "Missing: $file"
    echo "$file" >> "$MISSING_LOG"
    continue
  fi

  echo "Found: $file"

  # ROOT command (now includes 4th arg: IS_GREG_AI)
  cmd="clas12root -q 'source/gen_root_files_data.C(\"$hipopath\", \"$file\", \"${rootpath}${rootprefix}\", ${IS_GREG_AI})'"

  # Safe, readable session name + log file
  safe_file=$(basename "$file" | tr -c 'A-Za-z0-9_-' '_')
  sess="${JOB_PREFIX}_${safe_file}"
  LOG_FILE="${LOG_DIR}/${sess}.log"

  # Launch detached screen: export LD_LIBRARY_PATH, then run command, logging stdout+stderr
  screen -dmS "$sess" bash -lc "$LD_EXPORT; $cmd > \"$LOG_FILE\" 2>&1"

  printf "%s\t%s\t%s\t%s\n" "$sess" "$file" "$cmd" "$LOG_FILE" >> "$LAUNCHED_LOG"

  # Throttle to MAX_PARALLEL
  while (( $(running_jobs) >= MAX_PARALLEL )); do
    sleep 2
  done
done

echo
echo "Launched screen sessions with prefix: $JOB_PREFIX"
echo "Missing files logged to: $MISSING_LOG"
echo "Launched jobs logged to: $LAUNCHED_LOG"
echo "Per-session logs are in: $LOG_DIR"

cat <<EOF

Handy commands:
  screen -ls | grep $JOB_PREFIX
  screen -r <session_name>
  screen -S <session_name> -X quit
  tail -f logs/${JOB_PREFIX}_<sanitized-filename>.log

EOF