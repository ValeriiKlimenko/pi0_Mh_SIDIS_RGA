#!/usr/bin/env python3
"""
Run define_bin_migr_hists via ROOT CLI (link first, then execute).

You can either:
  * choose MX cut mode (recommended):
      --mx-cut-mode 0  -> rec_noMxcut/
      --mx-cut-mode 1  -> rec_Mxcut_1/
      --mx-cut-mode 2  -> rec_true/
  * or explicitly pass a directory with --dir (overrides mx-cut-mode)

Examples:
  python run_define_bin_migr_hists_cli.py
  python run_define_bin_migr_hists_cli.py --mx-cut-mode 0
  python run_define_bin_migr_hists_cli.py --mx-cut-mode 2
  python run_define_bin_migr_hists_cli.py --dir rec_noMxcut/
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

# Map MX cut modes to directory names
MX_DIR_MAP = {
    0: "rec_noMxcut/",
    1: "rec_Mxcut_1/",
    2: "rec_true/",
    3: "rec_match/",
}



def _q(s: str) -> str:
    """Escape for embedding in ROOT command lines."""
    return s.replace("\\", "\\\\").replace('"', '\\"')

def run_with_root_cli_link_first(dir_, outdir):
    if shutil.which("root") is None:
        raise SystemExit("ERROR: 'root' executable not found in PATH.")
    plus = "+" if COMPILE_MACRO else ""
    cmds = [
        f'.L {MACRO_PATH}{plus}',
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

    p.add_argument(
        "--mx-cut-mode",
        type=int,
        choices=[0, 1, 2],
        default=1,
        help=(
            "Missing mass cut mode: "
            "0 = no cut (rec_noMxcut/), "
            "1 = cut at 1.0 (rec_Mxcut_1/), "
            "2 = cut at 1.5 (rec_true/). "
            "Ignored if --dir is given."
        ),
    )

    p.add_argument(
        "--dir",
        default=None,
        help=(
            "Optional: explicit subdirectory under base data dir "
            "(e.g. rec_true/, rec_noMxcut/, rec_Mxcut_1/). "
            "If omitted, directory is chosen from --mx-cut-mode."
        ),
    )

    p.add_argument(
        "--isMCParentCut",
        type=int,
        default=0,
        help="Applies Parent/daugh. cut on MC rec gammas.",
    )

    # not implemented there
    p.add_argument(
        "--noPhiBinning",
        type=int,
        choices=[0, 1],
        default=0,
        help="Integrate over phi.",
    )

    args = p.parse_args()
    is_mc_parent_cuts = int(args.isMCParentCut)
    noPhiBinning = int(args.noPhiBinning)


    # Decide directory: explicit --dir wins, otherwise use mx-cut-mode mapping
    if args.dir is not None:
        dir_name = (args.dir or "").strip().strip("/") + "/"
    else:
        if (not is_mc_parent_cuts):
          dir_name = MX_DIR_MAP[args.mx_cut_mode]
        else:
          dir_name = MX_DIR_MAP[3]

    full_dir = os.path.join(BASE_DATA_DIR, dir_name)
    outdir = f"unfolding_{os.path.basename(os.path.normpath(dir_name))}"

    if not os.path.isdir(os.path.dirname(full_dir)) and not os.path.isdir(full_dir):
        print(f"Warning: constructed data dir does not exist: {full_dir}", file=sys.stderr)

    print(f"Base data dir : {BASE_DATA_DIR}")
    print(f"mx_cut_mode   : {args.mx_cut_mode}  -> {dir_name}")
    print(f"--dir         : {dir_name}")
    print(f"Full input dir: {full_dir}")
    print(f"Output folder : {outdir}")
    print(f"Macro         : {MACRO_PATH}  ({'compile +' if COMPILE_MACRO else 'interpret'})")

    run_with_root_cli_link_first(full_dir, outdir)
    print("Done.")

if __name__ == "__main__":
    main()
