#!/usr/bin/env python3
"""
Run define_bin_migr_hists via ROOT CLI (link first, then execute).
Only argument: --dir (rec_true/ or rec_fake/)

Examples:
  python run_define_bin_migr_hists_cli.py
  python run_define_bin_migr_hists_cli.py --dir rec_fake/
"""

import argparse
import os
import shutil
import subprocess
import sys

# Fixed settings
BASE_DATA_DIR = "/lustre24/expphy/volatile/clas12/valerii/multi_pi0/data_rec/"
MACRO_PATH    = "source/define_bin_migr_hists.cxx"
COMPILE_MACRO = False  # use ACLiC (+). Set False to interpret.

def _q(s: str) -> str:
    """Escape for embedding in ROOT command lines."""
    return s.replace("\\", "\\\\").replace('"', '\\"')

def run_with_root_cli_link_first(dir_, outdir):
    if shutil.which("root") is None:
        raise SystemExit("ERROR: 'root' executable not found in PATH.")
    plus = "+" if COMPILE_MACRO else ""
    cmds = [
        f'.L {MACRO_PATH}',
        f'define_bin_migr_hists("{_q(dir_)}","{_q(outdir)}")',
        '.q',
    ]
    cmd_stream = "\n".join(cmds) + "\n"

    print("ROOT command stream:\n" + "\n".join("  " + c for c in cmds))

    proc = subprocess.run(
        ["root", "-l", "-b"],
        input=cmd_stream,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: ROOT returned non-zero exit status {proc.returncode}")

def main():
    p = argparse.ArgumentParser(description="Run define_bin_migr_hists (ROOT CLI; link first).")
    p.add_argument("--dir", default="rec_true/", help="rec_true/ or rec_fake/ (relative under base data dir)")
    args = p.parse_args()

    # Normalize dir argument to always end with '/'
    dir_name = (args.dir or "").strip().strip("/") + "/"
    full_dir = os.path.join(BASE_DATA_DIR, dir_name)
    outdir = f"unfolding_{dir_name}"  # e.g. unfolding_rec_true/

    if not os.path.isdir(os.path.dirname(full_dir)) and not os.path.isdir(full_dir):
        print(f"Warning: constructed data dir does not exist: {full_dir}", file=sys.stderr)

    print(f"Base data dir : {BASE_DATA_DIR}")
    print(f"--dir         : {dir_name}")
    print(f"Full input dir: {full_dir}")
    print(f"Output folder : {outdir}")
    print(f"Macro         : {MACRO_PATH}  ({'compile +' if COMPILE_MACRO else 'interpret'})")

    run_with_root_cli_link_first(full_dir, outdir)
    print("Done.")

if __name__ == "__main__":
    main()
