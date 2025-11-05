#!/usr/bin/env python3
import argparse
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ----------------- CONFIG: set your prepared data directories here -----------------
DATA_DIR_MAP = {
    "DATA": "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_data/",
    "REC" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/",
    "GEN" : "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_gen/",
}
# ----------------------------------------------------------------------------------

MAX_SCREENS = 35
POLL_SEC = 1
SCREEN_NAME = "addBinning"

# Only queue data files (case-insensitive)
ALLOWED_EXTS = {".root", ".hipo"}  # tweak as needed

LOG = logging.getLogger("add_binning_launcher")

def setup_logging(logs_dir: Path, verbose: bool = True):
    logs_dir.mkdir(parents=True, exist_ok=True)
    master_log = logs_dir / "launcher.log"

    LOG.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(fmt)
    LOG.addHandler(ch)

    fh = logging.FileHandler(master_log, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    LOG.addHandler(fh)

    LOG.debug("Logging initialized. Master log: %s", master_log)

def run(cmd, *, cwd=None, check=True) -> subprocess.CompletedProcess:
    LOG.debug("Running: %s (cwd=%s)", cmd, cwd or os.getcwd())
    return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=check)

def screen_raw_ls() -> str:
    try:
        out = subprocess.check_output(["screen", "-ls"], text=True, stderr=subprocess.STDOUT)
        return out
    except subprocess.CalledProcessError as e:
        return e.output

def screen_ls():
    """
    Parse `screen -ls` and return {name: {pid, status}} for sessions that contain SCREEN_NAME.
    """
    out = screen_raw_ls()
    sessions = {}
    for line in out.splitlines():
        # lines look like: "    12345.job7   (Detached)"
        m = re.search(r"\s*(\d+)\.([^\s]+)\s+\(([^)]+)\)", line)
        if m:
            pid, name, status = int(m.group(1)), m.group(2), m.group(3)
            if SCREEN_NAME in name:
                sessions[name] = {"pid": pid, "status": status}
    return sessions

def verify_screen_started(name: str, timeout_sec: float = 10.0) -> bool:
    t0 = time.time()
    while time.time() - t0 < timeout_sec:
        if name in screen_ls():
            return True
        time.sleep(0.2)
    LOG.warning("Screen '%s' did not appear in `screen -ls` within %.1fs. Raw ls:\n%s",
                name, timeout_sec, screen_raw_ls())
    return False

def preflight_checks(logs_dir: Path):
    # Are tools present?
    for tool in ("screen", "bash"):
        if not shutil.which(tool):
            LOG.error("Required tool not found on PATH: %s", tool)
            raise SystemExit(1)
    if not shutil.which("root"):
        LOG.warning("`root` not found on PATH at launcher time. (It might still be available in the job environment.)")

    # Are we inside a screen already?
    if os.environ.get("STY"):
        LOG.info("You are currently inside a screen session (STY=%s). That's fine; just FYI.", os.environ["STY"])

    # Is logs dir writable?
    try:
        logs_dir.mkdir(parents=True, exist_ok=True)
        test_file = logs_dir / ".write_test"
        test_file.write_text("ok")
        test_file.unlink(missing_ok=True)
    except Exception as e:
        LOG.error("Logs directory is not writable: %s (%s)", logs_dir, e)
        raise SystemExit(1)

    # Macro exists?
    if not Path("source/add_binning.cxx").exists():
        LOG.error("Missing macro: source/add_binning.cxx")
        raise SystemExit(1)

def start_screen_session(name: str, bash_script: str, log_path: Path = None, workdir: Path = None):
    """
    Launch a detached screen session running the given bash_script.
    Output is appended to log_path if provided. Uses a finish() trap to write .OK/.FAIL.
    """
    # Header with strict mode + finish() trap
    header = f"""
set -Eeuo pipefail

_log={shlex.quote(str(log_path)) if log_path else ""}

finish() {{
  status=${{1:-$?}}
  # prevent double-run if both ERR and EXIT fire
  trap - ERR EXIT
  echo "[{name}] DONE at $(date) (exit=$status)"
  if [ -n "$_log" ]; then
    d=$(dirname "$_log")
    base=$(basename "$_log" .log)
    if [ $status -eq 0 ]; then
      touch "$d/$base.OK"
    else
      touch "$d/$base.FAIL"
    fi
  fi
  exit $status
}}

trap 'finish $?' ERR
trap 'finish $?' EXIT

echo "[{name}] START at $(date) (PID=$$)"
echo "[{name}] Working dir: $(pwd)"
echo "[{name}] STY=${{STY-}} (inside screen? $([ -n "${{STY-}}" ] && echo yes || echo no))"
"""

    # Build full script (no separate footer—finish() handles it)
    full_script = f"{header}\n{bash_script}\n"

    # One redirection at the end
    if log_path:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        redir = f">> {shlex.quote(str(log_path))} 2>&1"
    else:
        redir = ">/dev/null 2>&1"

    command = f"{full_script} {redir}"

    try:
        subprocess.run(
            ["screen", "-dmS", name, "bash", "-lc", command],
            cwd=str(workdir) if workdir else None,
            check=True,
            text=True,
        )
        LOG.info("Launched screen: %s (log: %s)", name, log_path if log_path else "<none>")
    except subprocess.CalledProcessError as e:
        LOG.error("Failed launching screen '%s': %s\nSTDOUT:\n%s\nSTDERR:\n%s",
                  name, e, e.stdout, e.stderr)
        raise

    if not verify_screen_started(name):
        LOG.error("Screen '%s' failed to start (not listed).", name)

if __name__ == "__main__":
    # --- CLI args ---
    parser = argparse.ArgumentParser(
        description="Launch add_binning ROOT macro in parallel screen sessions."
    )
    parser.add_argument(
        "--data-type",
        choices=["Data", "Rec", "Gen", "DATA", "REC", "GEN"],
        default="Data",
        help="Which dataset to process: Data, Rec, or Gen (case-insensitive).",
    )
    parser.add_argument(
        "--is-true-gen-event",
        type=int,
        choices=[0, 1],
        default=1,
        help="Pass 0 or 1 to the macro as the is_true_gen_event flag (default: 1).",
    )
    parser.add_argument(
        "--logs-dir",
        default="logs_add_binning",
        help="Directory to store screen logs (also stores per-job .OK/.FAIL sentinels and launcher.log).",
    )
    parser.add_argument(
        "--bin-test",
        type=int,
        choices=[0, 1],
        default=0,
        help="fills z, pt2, phi in the output ttree (0/1).",
    )
    parser.add_argument(
        "--max-screens",
        type=int,
        default=MAX_SCREENS,
        help="Maximum concurrent screens.",
    )
    parser.add_argument(
        "--poll-sec",
        type=float,
        default=POLL_SEC,
        help="Polling interval when waiting for free screen slots.",
    )
    args = parser.parse_args()

    logs_dir = Path(args.logs_dir)
    setup_logging(logs_dir)
    preflight_checks(logs_dir)

    # Normalize data type and pick directory/tag
    data_tag = args.data_type.upper()  # "DATA", "REC", "GEN"
    if data_tag not in DATA_DIR_MAP:
        LOG.error("Unknown data type '%s'. Expected one of: Data, Rec, Gen.", args.data_type)
        raise SystemExit(1)
    data_dir = DATA_DIR_MAP[data_tag]
    LOG.info("Selected data type: %s  -> %s", data_tag, data_dir)

    # Validate input dir
    mypath = Path(data_dir)
    if not mypath.exists() or not mypath.is_dir():
        LOG.error("Data directory does not exist or is not a directory: %s", mypath)
        raise SystemExit(1)

    # Gather input files (case-insensitive by suffix)
    rec_files = sorted(
        [p.name for p in mypath.iterdir() if p.is_file() and p.suffix.lower() in ALLOWED_EXTS]
    )
    if not rec_files:
        LOG.error("No input data files (%s) found in: %s", ",".join(ALLOWED_EXTS), mypath)
        raise SystemExit(1)

    LOG.debug("First few files: %s", rec_files[:5])

    is_true_gen_event = int(args.is_true_gen_event)
    is_bin_test = int(args.bin_test)
    LOG.info("Files to process: %d  | is_true_gen_event=%d  | bin_test=%d",
             len(rec_files), is_true_gen_event, is_bin_test)

    # ---- UPDATED LIB PATHS ----
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

    try:
        for i_file, file in enumerate(rec_files):
            # wait until there is an available screen slot
            while len(screen_ls()) >= args.max_screens:
                LOG.debug("Max screens reached (%d). Sleeping %.2fs ...", args.max_screens, args.poll_sec)
                time.sleep(args.poll_sec)

            file_path = str(mypath / file)
            name = f"{SCREEN_NAME}_{i_file}"           # unique screen session name
            log = logs_dir / f"{name}.log"             # capture stdout/stderr

            # ROOT macro call
            macro_call = f'source/add_binning.cxx("{file_path}", "{data_tag}", {is_true_gen_event}, {is_bin_test})'
            root_call = f"root -l -b -q {shlex.quote(macro_call)}"

            # Per-job bash script
            bash_script = f"""
echo "[{name}] Host: $(hostname)"
echo "[{name}] Start: $(date)"
{LD_EXPORT}
umask 002
echo "[{name}] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "[{name}] Running: {root_call}"

set +e
{root_call}
ec=$?
set -e

echo "[{name}] ROOT exit code: $ec"
exit $ec
"""

            start_screen_session(name, bash_script, log_path=log)
            LOG.info("Launched %s [data_type=%s] (is_true_gen_event=%d, bin_test=%d)  (log: %s)",
                     name, data_tag, is_true_gen_event, is_bin_test, log)

    except KeyboardInterrupt:
        LOG.warning("Interrupted by user. Current screens:\n%s", screen_raw_ls())
        raise
