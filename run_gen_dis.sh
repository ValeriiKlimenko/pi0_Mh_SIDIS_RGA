#!/usr/bin/env bash
# Run gen_root_files_rec.C jobs in parallel inside detached GNU screen sessions
# with LD_LIBRARY_PATH set per session and logs per job.

# --- user settings ---
input="path_to_files/nSidis_45_50nA_092025_v2.dat"
hipopath="/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/"
rootpath="/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen_dis/"
rootprefix="gen_dis_f2018_"

# How many screen sessions to run at once
MAX_PARALLEL=20

# --- safety checks & prep ---
command -v screen >/dev/null 2>&1 || { echo "ERROR: 'screen' not found in PATH." >&2; exit 1; }
[[ -f "$input" ]] || { echo "ERROR: Input file list not found: $input" >&2; exit 1; }

JOB_PREFIX="gen_dis_$(date +%Y%m%d_%H%M%S)"
MISSING_LOG="gen_dis_file.log"
LAUNCHED_LOG="launched_gen_dis_jobs.log"
LOG_DIR="logs"

mkdir -p "$LOG_DIR"
: > "$MISSING_LOG"
: > "$LAUNCHED_LOG"

# Load & normalize file list (strip hipopath prefix if present)
mapfile -t raw_lines < "$input"
files=()
for line in "${raw_lines[@]}"; do
  [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
  line="${line#$hipopath}"
  files+=("$line")
done

echo "Read ${#files[@]} gen dis (clasdis) files"

# Helper: count running screen jobs for this batch
running_jobs() {
  screen -ls 2>/dev/null | grep -c "$JOB_PREFIX" || true
}

# LD_LIBRARY_PATH export to run inside each screen
LD_EXPORT='export LD_LIBRARY_PATH="/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/root/6.30.04/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/clas12root/1.8.5/4.2.0/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/ccdb/1.99.6/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/iguana/0.8.0/4.2.0/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0/lib:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib64:/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/lib:${LD_LIBRARY_PATH-}"'

# --- main loop ---
for file in "${files[@]}"; do
  fpath="${hipopath}${file}"
  if [[ ! -f "$fpath" ]]; then
    echo "Missing: $file"
    echo "$file" >> "$MISSING_LOG"
    continue
  fi

  echo "Found: $file"

  # ROOT command
  cmd="clas12root -q 'source/gen_root_files_gen_dis.C(\"$hipopath\", \"$file\", \"${rootpath}${rootprefix}\")'"

  # Screen session name + logfile
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

